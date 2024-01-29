"""
calculate_SharedMutation_gene.py: calculate the number of of shared mutations in a gene level
"""
import argparse
import multiprocessing
import tarfile
import typing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import tqdm
import step00

lower_precancer_list: typing.List[str] = list()
higher_precancer_list: typing.List[str] = list()

lower_patients_length = 0
higher_patients_length = 0

data = pandas.DataFrame()


def query_mutation(gene: str, sample: str) -> int:
    return len(data.loc[(data["Precancer"] == sample) & (data["Hugo_Symbol"] == gene)])


def query_exact(gene: str, derivation: str) -> typing.Optional[float]:
    confusion_matrix = numpy.zeros((2, 2))

    try:
        confusion_matrix[0, 0] = len(set(data.loc[(data["Precancer"].isin(lower_precancer_list)) & (data["Hugo_Symbol"] == gene), "Patient"]))
    except KeyError:
        pass

    confusion_matrix[0, 1] = lower_patients_length - confusion_matrix[0, 0]

    try:
        confusion_matrix[1, 0] = len(set(data.loc[(data["Precancer"].isin(higher_precancer_list)) & (data["Hugo_Symbol"] == gene), "Patient"]))
    except KeyError:
        pass

    confusion_matrix[1, 1] = higher_patients_length - confusion_matrix[1, 0]

    if derivation == "Fisher":
        return scipy.stats.fisher_exact(confusion_matrix)[1]
    elif derivation == "Barnard":
        return scipy.stats.barnard_exact(confusion_matrix).pvalue
    elif derivation == "Boschloo":
        return scipy.stats.boschloo_exact(confusion_matrix).pvalue
    else:
        raise Exception("Something went wrong!!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Shared mutation information TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data w/ Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--percentage", help="Percentage of patients to include", type=float, default=0.1)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.05)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0.0 < args.percentage < 0.5):
        raise ValueError("Percentage must be (0.0, 0.5)!!")
    elif not (0 < args.p < 1):
        raise ValueError("P-values must be (0, 1)")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(len(patients), sorted(patients))

    input_data = input_data.loc[(input_data["Patient"].isin(patients))]
    print(input_data)

    filtered_data = input_data.loc[(input_data[step00.nonsynonymous_column].isin(step00.nonsynonymous_mutations))]
    print(filtered_data)

    gene_list = sorted(set(input_data["Hugo_Symbol"]))
    print("Gene:", len(gene_list))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    for MSP in tqdm.tqdm(step00.sharing_columns):
        clinical_data = clinical_data.sort_values(MSP)

        lower_bound, higher_bound = numpy.quantile(clinical_data[MSP], args.percentage), numpy.quantile(clinical_data[MSP], 1.0 - args.percentage)

        lower_precancer_list = list(clinical_data.loc[(clinical_data[MSP] < lower_bound), f"{MSP}-sample"])
        higher_precancer_list = list(clinical_data.loc[(clinical_data[MSP] > higher_bound), f"{MSP}-sample"])

        lower_patients = list(map(step00.get_patient, lower_precancer_list))
        higher_patients = list(map(step00.get_patient, higher_precancer_list))

        lower_patients_length = len(lower_patients)
        higher_patients_length = len(higher_patients)

        if "SYN" in MSP:
            data = input_data.loc[(input_data["Patient"].isin(lower_patients + higher_patients))]
        else:
            data = filtered_data.loc[(filtered_data["Patient"].isin(lower_patients + higher_patients))]

        heatmap_data = pandas.DataFrame(data=numpy.zeros((len(gene_list), len(lower_precancer_list + higher_precancer_list))), index=gene_list, columns=lower_precancer_list + higher_precancer_list, dtype=int)
        with multiprocessing.Pool(args.cpus) as pool:
            for sample in tqdm.tqdm(list(heatmap_data.columns), leave=False):
                heatmap_data[sample] = list(pool.starmap(query_mutation, [(gene, sample) for gene in gene_list]))

        exact_test_data = pandas.DataFrame(data=numpy.zeros((len(gene_list), 3)), index=gene_list, columns=["Fisher", "Barnard", "Boschloo"], dtype=float)
        with multiprocessing.Pool(args.cpus) as pool:
            for derivation in tqdm.tqdm(list(exact_test_data.columns), leave=False):
                exact_test_data[derivation] = -1 * numpy.log10(numpy.array(pool.starmap(query_exact, [(gene, derivation) for gene in list(exact_test_data.index)])))

        exact_test_data = exact_test_data.loc[(exact_test_data > (-1 * numpy.log10(args.p))).all(axis="columns")].sort_values(by="Fisher", kind="stable", ascending=False)
        heatmap_data = heatmap_data.loc[exact_test_data.index, :]

        fig, axs = matplotlib.pyplot.subplots(ncols=3, figsize=(24, 24), gridspec_kw={"width_ratios": [lower_patients_length + 5, len(exact_test_data.columns) + 5, higher_patients_length + 5]})

        seaborn.heatmap(data=heatmap_data.loc[:, lower_precancer_list], vmin=0, vmax=heatmap_data.max().max(), cmap="gray", cbar=False, xticklabels=lower_patients, yticklabels=True, fmt="d", annot=True, ax=axs[0])
        axs[0].set_xlabel("Lower")

        seaborn.heatmap(data=exact_test_data, cmap="Reds", vmin=0, center=-1 * numpy.log10(args.p), cbar=True, xticklabels=True, yticklabels=True, ax=axs[1])

        seaborn.heatmap(data=heatmap_data.loc[:, higher_precancer_list], vmin=0, vmax=heatmap_data.max().max(), cmap="gray", cbar=False, xticklabels=higher_patients, yticklabels=True, fmt="d", annot=True, ax=axs[2])
        axs[2].set_xlabel("Higher")

        matplotlib.pyplot.tight_layout()
        figures.append(f"{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(figures):
            tar.add(f, arcname=f)
