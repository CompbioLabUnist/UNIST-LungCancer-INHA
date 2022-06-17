"""
plot_gene_sharing.py: Plot exact test for gene in mutation sharing
"""
import argparse
import collections
import multiprocessing
import typing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import tqdm
import step00

mutect_data = pandas.DataFrame()
mutation_set: collections.Counter = collections.Counter()
heatmap_data = pandas.DataFrame()


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


def query_mutect(gene: str, sample: str) -> int:
    return mutation_set[(gene, sample)]


def query_heatmap(gene: str, derivation: str) -> typing.Union[float, None]:
    confusion_matrix = numpy.zeros((2, 2))
    try:
        confusion_matrix[0, 0] = heatmap_data.loc[gene, lower_patients].astype(bool).T.value_counts().loc[True]
    except KeyError:
        pass

    try:
        confusion_matrix[0, 1] = heatmap_data.loc[gene, lower_patients].astype(bool).T.value_counts().loc[False]
    except KeyError:
        pass

    try:
        confusion_matrix[1, 0] = heatmap_data.loc[gene, higher_patients].astype(bool).T.value_counts().loc[True]
    except KeyError:
        pass

    try:
        confusion_matrix[1, 1] = heatmap_data.loc[gene, higher_patients].astype(bool).T.value_counts().loc[False]
    except KeyError:
        pass

    if derivation == "Fisher":
        return scipy.stats.fisher_exact(confusion_matrix)[1]
    elif derivation == "Chi2":
        return scipy.stats.chi2_contingency(confusion_matrix)[1]
    elif derivation == "Barnard":
        return scipy.stats.barnard_exact(confusion_matrix).pvalue
    elif derivation == "Boschloo":
        return scipy.stats.boschloo_exact(confusion_matrix).pvalue
    else:
        raise Exception("Something went wrong!!")


def query_mutation(gene: str, patient: str) -> int:
    wanted_columns = ["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]

    patient_data = mutect_data.loc[(mutect_data["Patient"] == patient) & (mutect_data["Hugo_Symbol"] == gene)]

    stage_set = list(filter(lambda x: x in set(patient_data["Stage"]), step00.long_sample_type_list))
    if ("Primary" in stage_set) and (len(stage_set) > 1):
        primary_set = set(patient_data.loc[patient_data["Stage"] == "Primary", wanted_columns].itertuples(index=False, name=None))
        precancer_set = set(patient_data.loc[patient_data["Stage"] == stage_set[-2], wanted_columns].itertuples(index=False, name=None))
        return len(precancer_set & primary_set)
    else:
        return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("cgc", help="CGC gene CSV file", type=str)
    parser.add_argument("tsv", help="Output TSV file", type=str)
    parser.add_argument("pdf", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--threshold", help="Threshold to use", type=int, default=100)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.05)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_strategy = parser.add_mutually_exclusive_group(required=True)
    group_strategy.add_argument("--median", help="Median division", action="store_true", default=False)
    group_strategy.add_argument("--mean", help="Mean division", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.cgc.endswith(".csv"):
        raise ValueError("CGC must end with .CSV!!")
    elif not args.tsv.endswith(".tsv"):
        raise ValueError("TSV must end with .tsv!!")
    elif not args.pdf.endswith(".pdf"):
        raise ValueError("PDF must end with .pdf!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.p < 1):
        raise ValueError("P-values must be (0, 1)")

    clinical_data: pandas.DataFrame = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(patients)

    cgc_data = pandas.read_csv(args.cgc, index_col="Gene Symbol")
    print(cgc_data)

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))
    args.input.sort(key=step00.sorting)
    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
        mutect_data["Tumor_Sample_Barcode"] = pool.map(step00.get_id, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Patient"] = pool.map(step00.get_patient, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Stage"] = pool.map(step00.get_long_sample_type, mutect_data["Tumor_Sample_Barcode"])
    print(mutect_data)

    mutect_data = mutect_data[(mutect_data["Variant_Classification"].isin(step00.nonsynonymous_mutations))]
    print(mutect_data)

    patients &= set(mutect_data["Patient"])

    wanted_columns = ["Hugo_Symbol", "Patient", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]
    clinical_data["Shared Proportion"] = None
    for patient in tqdm.tqdm(patients):
        patient_data = mutect_data.loc[mutect_data["Patient"] == patient]

        stage_set = list(filter(lambda x: x in set(patient_data["Stage"]), step00.long_sample_type_list))
        assert "Primary" in stage_set
        primary_set = set(patient_data.loc[patient_data["Stage"] == "Primary", wanted_columns].itertuples(index=False, name=None))

        precancer_set = set(patient_data.loc[patient_data["Stage"] == stage_set[-2], wanted_columns].itertuples(index=False, name=None))
        clinical_data.loc[patient, "Shared Proportion"] = len(primary_set & precancer_set) / len(primary_set)
        mutation_set += collections.Counter(list(map(lambda x: x[:2], primary_set & precancer_set)))
    clinical_data.dropna(subset=["Shared Proportion"], inplace=True)
    print(mutation_set.most_common(5))
    print(clinical_data)

    if args.median:
        cutting = numpy.median(clinical_data["Shared Proportion"])
    elif args.mean:
        cutting = numpy.mean(clinical_data["Shared Proportion"])
    else:
        raise Exception("Something went wrong!!")

    lower_data = clinical_data.loc[(clinical_data["Shared Proportion"] <= cutting)].sort_values(by="Shared Proportion")
    lower_patients = list(lower_data.index)
    higher_data = clinical_data.loc[(clinical_data["Shared Proportion"] > cutting)].sort_values(by="Shared Proportion")
    higher_patients = list(higher_data.index)
    print(len(lower_patients), "vs", len(higher_patients))

    gene_list = sorted(set(cgc_data.index) & set(map(lambda x: x[0], mutation_set.keys())))
    print(len(gene_list))

    heatmap_data = pandas.DataFrame(data=numpy.zeros((len(gene_list), len(lower_patients) + len(higher_patients))), index=gene_list, columns=lower_patients + higher_patients, dtype=int)
    with multiprocessing.Pool(args.cpus) as pool:
        for sample in tqdm.tqdm(list(heatmap_data.columns)):
            heatmap_data.loc[:, sample] = pool.starmap(query_mutect, [(gene, sample) for gene in gene_list])
    print(heatmap_data)

    mutation_data = pandas.DataFrame(index=gene_list, columns=lower_patients + higher_patients, dtype=int)
    with multiprocessing.Pool(args.cpus) as pool:
        for sample in tqdm.tqdm(list(mutation_data.columns)):
            mutation_data.loc[:, sample] = pool.starmap(query_mutation, [(gene, sample) for gene in gene_list])
    print(mutation_data)

    exact_test_data = pandas.DataFrame(data=numpy.zeros((len(gene_list), 4)), index=gene_list, columns=["Fisher", "Chi2", "Barnard", "Boschloo"], dtype=float)
    with multiprocessing.Pool(args.cpus) as pool:
        for derivation in tqdm.tqdm(list(exact_test_data.columns)):
            exact_test_data.loc[:, derivation] = -1 * numpy.log10(numpy.array(pool.starmap(query_heatmap, [(gene, derivation) for gene in list(exact_test_data.index)])))
    print(exact_test_data)

    exact_test_data = exact_test_data.loc[(exact_test_data > -1 * numpy.log10(args.p)).any(axis="columns")].sort_values(by="Fisher", kind="stable", ascending=False)
    heatmap_data = heatmap_data.loc[exact_test_data.index, :]
    mutation_data = mutation_data.loc[exact_test_data.index, :]
    print(exact_test_data)

    mutation_data.columns = list(map(lambda x: f"{x}-Lower" if (x in lower_patients) else f"{x}-Higher", list(mutation_data.columns)))
    output_data = pandas.concat([exact_test_data, mutation_data], axis="columns", join="outer", verify_integrity=True)
    output_data.to_csv(args.tsv, sep="\t")
    print(output_data)

    exact_test_data = exact_test_data.iloc[:args.threshold, :]
    heatmap_data = heatmap_data.loc[exact_test_data.index, :]
    mutation_data = mutation_data.loc[exact_test_data.index, :]

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, axs = matplotlib.pyplot.subplots(ncols=3, figsize=(len(lower_patients) + len(exact_test_data.columns) + len(higher_patients) + 15, exact_test_data.shape[0] + 5), gridspec_kw={"width_ratios": [len(lower_patients) + 5, len(exact_test_data.columns) + 5, len(higher_patients) + 5]})

    seaborn.heatmap(data=heatmap_data.loc[:, lower_patients], vmin=0, vmax=heatmap_data.max().max(), cmap="gray", cbar=False, xticklabels=True, yticklabels=True, fmt="d", annot=True, ax=axs[0])
    axs[0].set_xlabel(f"Lower {len(lower_patients)} Patients")

    seaborn.heatmap(data=exact_test_data, cmap="Reds", vmin=0, center=-1 * numpy.log10(args.p), cbar=True, xticklabels=True, yticklabels=True, ax=axs[1])

    seaborn.heatmap(data=heatmap_data.loc[:, higher_patients], vmin=0, vmax=heatmap_data.max().max(), cmap="gray", cbar=False, xticklabels=True, yticklabels=True, fmt="d", annot=True, ax=axs[2])
    axs[2].set_xlabel(f"Higher {len(higher_patients)} Patients")

    matplotlib.pyplot.tight_layout()
    fig.savefig(args.pdf)
    matplotlib.pyplot.close(fig)
