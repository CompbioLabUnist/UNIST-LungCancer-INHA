"""
plot_gene_clinical.py: Plot the importance of gene upon clinical data
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
        confusion_matrix[0, 0] = heatmap_data.loc[gene, control_samples].astype(bool).T.value_counts().loc[True]
    except KeyError:
        pass

    try:
        confusion_matrix[0, 1] = heatmap_data.loc[gene, control_samples].astype(bool).T.value_counts().loc[False]
    except KeyError:
        pass

    try:
        confusion_matrix[1, 0] = heatmap_data.loc[gene, case_samples].astype(bool).T.value_counts().loc[True]
    except KeyError:
        pass

    try:
        confusion_matrix[1, 1] = heatmap_data.loc[gene, case_samples].astype(bool).T.value_counts().loc[False]
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


def query_mutation(gene: str, sample: str) -> str:
    return ",".join(sorted(mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == sample) & (mutect_data["Hugo_Symbol"] == gene), "Variant_Classification"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("table", help="Output TSV file", type=str)
    parser.add_argument("figure", help="Output PDF file", type=str)
    parser.add_argument("--compare", help="Comparison grouping (type, control, case)", type=str, nargs=3, default=["Recurrence", "NO", "YES"])
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.05)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.table.endswith(".tsv"):
        raise ValueError("Table must end with .TSV!!")
    elif not args.figure.endswith(".pdf"):
        raise ValueError("Figure must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.p < 1):
        raise ValueError("P-values must be (0, 1)")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        control_patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC") & (clinical_data[args.compare[0]] == args.compare[1])].index)
        case_patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC") & (clinical_data[args.compare[0]] == args.compare[2])].index)
    elif args.ADC:
        control_patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC") & (clinical_data[args.compare[0]] == args.compare[1])].index)
        case_patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC") & (clinical_data[args.compare[0]] == args.compare[2])].index)
    else:
        raise Exception("Something went wrong!!")
    print(sorted(control_patients))
    print(sorted(case_patients))

    args.input = sorted(filter(lambda x: step00.get_patient(x.split("/")[-1].split(".")[0]) in control_patients, args.input), key=step00.sorting_by_type) + sorted(filter(lambda x: step00.get_patient(x.split("/")[-1].split(".")[0]) in case_patients, args.input), key=step00.sorting_by_type)

    control_samples = list(map(step00.get_id, sorted(filter(lambda x: step00.get_patient(x.split("/")[-1].split(".")[0]) in control_patients, args.input), key=step00.sorting_by_type)))
    case_samples = list(map(step00.get_id, sorted(filter(lambda x: step00.get_patient(x.split("/")[-1].split(".")[0]) in case_patients, args.input), key=step00.sorting_by_type)))

    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)

    mutect_data = mutect_data.loc[(mutect_data["Variant_Classification"].isin(step00.nonsynonymous_mutations))]
    mutect_data["Tumor_Sample_Barcode"] = list(map(lambda x: x.split(".")[0], mutect_data["Tumor_Sample_Barcode"]))
    print(mutect_data)

    mutation_set = collections.Counter(mutect_data[["Hugo_Symbol", "Tumor_Sample_Barcode"]].itertuples(index=False, name=None))
    gene_list = sorted(set(mutect_data["Hugo_Symbol"]))

    heatmap_data = pandas.DataFrame(data=numpy.zeros((len(gene_list), len(args.input))), index=gene_list, columns=control_samples + case_samples, dtype=int)
    with multiprocessing.Pool(args.cpus) as pool:
        for gene in tqdm.tqdm(list(heatmap_data.index)):
            heatmap_data.loc[gene, :] = pool.starmap(query_mutect, [(gene, sample) for sample in (control_samples + case_samples)])
    print(heatmap_data)

    mutation_data = pandas.DataFrame(index=gene_list, columns=control_samples + case_samples, dtype=str)
    with multiprocessing.Pool(args.cpus) as pool:
        for sample in tqdm.tqdm(list(heatmap_data.columns)):
            heatmap_data.loc[:, sample] = pool.starmap(query_mutect, [(gene, sample) for gene in gene_list])
    print(mutation_data)

    exact_test_data = pandas.DataFrame(data=numpy.zeros((len(gene_list), 4)), index=gene_list, columns=["Fisher", "Chi2", "Barnard", "Boschloo"], dtype=float)
    with multiprocessing.Pool(args.cpus) as pool:
        for derivation in tqdm.tqdm(list(exact_test_data.columns)):
            exact_test_data.loc[:, derivation] = -1 * numpy.log10(pool.starmap(query_heatmap, [(gene, derivation) for gene in list(exact_test_data.index)]))
    print(exact_test_data)

    exact_test_data = exact_test_data.loc[(exact_test_data > -1 * numpy.log10(args.p)).any(axis="columns")].sort_values(by="Fisher", ascending=False).iloc[:100, :]
    print(exact_test_data)

    heatmap_data = heatmap_data.loc[exact_test_data.index, :]

    fig, axs = matplotlib.pyplot.subplots(ncols=3, figsize=(len(control_samples) + len(exact_test_data.columns) + len(case_samples) + 15, exact_test_data.shape[0] + 5), gridspec_kw={"width_ratios": [len(control_samples) + 5, len(exact_test_data.columns) + 5, len(case_samples) + 5]})

    seaborn.heatmap(data=heatmap_data.loc[:, control_samples], vmin=0, vmax=heatmap_data.max().max(), cmap="gray", cbar=False, xticklabels=True, yticklabels=True, fmt="d", annot=True, ax=axs[0])
    axs[0].set_xlabel("{0} - {1}".format(args.compare[0], args.compare[1]))

    seaborn.heatmap(data=exact_test_data, cmap="Reds", vmin=0, center=-1 * numpy.log10(args.p), cbar=True, xticklabels=True, yticklabels=True, ax=axs[1])

    seaborn.heatmap(data=heatmap_data.loc[:, case_samples], vmin=0, vmax=heatmap_data.max().max(), cmap="gray", cbar=False, xticklabels=True, yticklabels=True, fmt="d", annot=True, ax=axs[2])
    axs[2].set_xlabel("{0} - {1}".format(args.compare[0], args.compare[2]))

    matplotlib.pyplot.tight_layout()
    fig.savefig(args.figure)
    matplotlib.pyplot.close(fig)

    mutation_data = mutation_data.loc[:, sorted(mutation_data.columns, key=step00.sorting)]
    mutation_data.columns = list(map(lambda x: "{0}-{1}".format(x, args.compare[1]) if x in control_samples else "{0}-{1}".format(x, args.compare[2]), list(mutation_data.columns)))
    output_data = pandas.concat([exact_test_data, mutation_data], axis="columns", join="outer", verify_integrity=True)
    output_data.to_csv(args.table, sep="\t")
    print(output_data)
