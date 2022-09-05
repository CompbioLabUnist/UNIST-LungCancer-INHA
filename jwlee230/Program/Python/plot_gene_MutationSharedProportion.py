"""
plot_gene_MutationSharedProportion.py: plot gene exact test upon Mutation Shared Proportion
"""
import argparse
import collections
import itertools
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

mutect_data = pandas.DataFrame()
filtered_mutect_data = pandas.DataFrame()
data = pandas.DataFrame()
mutation_set: collections.Counter = collections.Counter()
heatmap_data = pandas.DataFrame()

lower_samples: typing.List[str] = list()
higher_samples: typing.List[str] = list()


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


def query_mutect(gene: str, sample: str) -> int:
    return mutation_set[(gene, sample)]


def query_heatmap(gene: str, derivation: str) -> typing.Union[float, None]:
    confusion_matrix = numpy.zeros((2, 2))
    try:
        confusion_matrix[0, 0] = heatmap_data.loc[gene, lower_samples].astype(bool).T.value_counts().loc[True]
    except KeyError:
        pass

    try:
        confusion_matrix[0, 1] = heatmap_data.loc[gene, lower_samples].astype(bool).T.value_counts().loc[False]
    except KeyError:
        pass

    try:
        confusion_matrix[1, 0] = heatmap_data.loc[gene, higher_samples].astype(bool).T.value_counts().loc[True]
    except KeyError:
        pass

    try:
        confusion_matrix[1, 1] = heatmap_data.loc[gene, higher_samples].astype(bool).T.value_counts().loc[False]
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
    return ",".join(sorted(list(map(lambda x: step00.nonsynonymous_notations[x], data.loc[(data["Tumor_Sample_Barcode"] == sample) & (data["Hugo_Symbol"] == gene), "Variant_Classification"]))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("cgc", help="CGC gene CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
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
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.cgc.endswith(".csv"):
        raise ValueError("CGC must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Figure must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif args.threshold <= 50:
        raise ValueError("Threshold must be greater than 50!!")
    elif not (0 < args.p < 1):
        raise ValueError("P-values must be (0, 1)")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    cgc_data = pandas.read_csv(args.cgc, index_col="Gene Symbol")
    cgc_data = cgc_data.loc[~(cgc_data["Genome Location"].str.contains(":-"))]
    cgc_data["chromosome"] = list(map(lambda x: "chr" + (x.replace("-", ":").split(":")[0]), cgc_data["Genome Location"]))
    cgc_data["start"] = list(map(lambda x: int(x.replace("-", ":").split(":")[1]), cgc_data["Genome Location"]))
    cgc_data["end"] = list(map(lambda x: int(x.replace("-", ":").split(":")[2]), cgc_data["Genome Location"]))
    print(cgc_data)

    gene_list = list(cgc_data.index)

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(patients)

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))
    args.input.sort(key=step00.sorting)
    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
        mutect_data["Tumor_Sample_Barcode"] = pool.map(step00.get_id, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Patient"] = pool.map(step00.get_patient, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Stage"] = pool.map(step00.get_long_sample_type, mutect_data["Tumor_Sample_Barcode"])
    print(mutect_data)

    mutect_data = mutect_data.loc[(mutect_data["Hugo_Symbol"].isin(gene_list))]
    print(mutect_data)

    for MSP in tqdm.tqdm(step00.sharing_columns):
        if args.median:
            cutting = numpy.median(clinical_data[MSP])
        elif args.mean:
            cutting = numpy.mean(clinical_data[MSP])
        else:
            raise Exception("Something went wrong!!")

        selecting = set(clinical_data.loc[(clinical_data[MSP] <= cutting)].index)
        mutect_data[MSP] = list(map(lambda x: "Lower" if (step00.get_patient(x) in selecting) else "Higher", mutect_data["Tumor_Sample_Barcode"]))
    print(mutect_data)

    filtered_mutect_data = mutect_data[(mutect_data[step00.nonsynonymous_column].isin(step00.nonsynonymous_mutations))]
    print(filtered_mutect_data)

    stage_list = list(filter(lambda x: x in set(mutect_data["Stage"]), step00.long_sample_type_list))
    print(stage_list)

    files: typing.List[str] = list()
    for MSP, stage in tqdm.tqdm(list(itertools.product(step00.sharing_columns[:1], stage_list))):
        if "SYN" in MSP:
            data = mutect_data.loc[(mutect_data["Stage"] == stage)]
        else:
            data = filtered_mutect_data.loc[(filtered_mutect_data["Stage"] == stage)]

        if args.median:
            cutting = numpy.median(clinical_data[MSP])
        elif args.mean:
            cutting = numpy.mean(clinical_data[MSP])
        else:
            raise Exception("Something went wrong!!")

        lower_samples = sorted(set(data.loc[(data[MSP] == "Lower"), "Tumor_Sample_Barcode"]), key=lambda x: clinical_data.loc[step00.get_patient(x), MSP])
        higher_samples = sorted(set(data.loc[(data[MSP] == "Higher"), "Tumor_Sample_Barcode"]), key=lambda x: clinical_data.loc[step00.get_patient(x), MSP])

        if (not lower_samples) or (not higher_samples):
            continue

        mutation_set = collections.Counter(data[["Hugo_Symbol", "Tumor_Sample_Barcode"]].itertuples(index=False, name=None))
        gene_list = sorted(set(data["Hugo_Symbol"]))

        heatmap_data = pandas.DataFrame(data=numpy.zeros((len(gene_list), len(lower_samples + higher_samples))), index=gene_list, columns=lower_samples + higher_samples, dtype=int)
        with multiprocessing.Pool(args.cpus) as pool:
            for sample in tqdm.tqdm(list(heatmap_data.columns)):
                heatmap_data.loc[:, sample] = pool.starmap(query_mutect, [(gene, sample) for gene in gene_list])
        print(heatmap_data)

        mutation_data = pandas.DataFrame(index=gene_list, columns=lower_samples + higher_samples, dtype=str)
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

        mutation_data.columns = list(map(lambda x: f"{x}-Lower" if (x in lower_samples) else f"{x}-Higher", list(mutation_data.columns)))
        output_data = pandas.concat([exact_test_data, mutation_data], axis="columns", join="outer", verify_integrity=True)
        files.append(f"{stage}_{MSP}.tsv")
        output_data.to_csv(files[-1], sep="\t")
        print(output_data)

        exact_test_data = exact_test_data.iloc[:args.threshold, :]
        heatmap_data = heatmap_data.loc[exact_test_data.index, :]
        mutation_data = mutation_data.loc[exact_test_data.index, :]

        fig, axs = matplotlib.pyplot.subplots(ncols=3, figsize=(len(lower_samples) + len(exact_test_data.columns) + len(higher_samples) + 15, exact_test_data.shape[0] + 5), gridspec_kw={"width_ratios": [len(lower_samples) + 5, len(exact_test_data.columns) + 5, len(higher_samples) + 5]})

        seaborn.heatmap(data=heatmap_data.loc[:, lower_samples], vmin=0, vmax=heatmap_data.max().max(), cmap="gray", cbar=False, xticklabels=True, yticklabels=True, fmt="d", annot=True, ax=axs[0])
        axs[0].set_xlabel(f"{MSP} - Lower")

        seaborn.heatmap(data=exact_test_data, cmap="Reds", vmin=0, center=-1 * numpy.log10(args.p), cbar=True, xticklabels=True, yticklabels=True, ax=axs[1])

        seaborn.heatmap(data=heatmap_data.loc[:, higher_samples], vmin=0, vmax=heatmap_data.max().max(), cmap="gray", cbar=False, xticklabels=True, yticklabels=True, fmt="d", annot=True, ax=axs[2])
        axs[2].set_xlabel(f"{MSP} - Higher")

        matplotlib.pyplot.tight_layout()
        files.append(f"{stage}_{MSP}.pdf")
        fig.savefig(files[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(files):
            tar.add(f, arcname=f)
