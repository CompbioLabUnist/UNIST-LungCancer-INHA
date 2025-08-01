"""
draw_vaf_plots.py: Draw VAF plots upon Mutect2 MAF results
"""
import argparse
import multiprocessing
import matplotlib
import matplotlib.pyplot
import seaborn
import sklearn.manifold
import sklearn.preprocessing
import pandas
import step00

gene_set = set()


def read(filename: str) -> pandas.DataFrame:
    sample = filename.split("/")[-1].split(".")[0]

    data = pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)

    # data = data.loc[(data["Chromosome"].isin(step00.chromosome_list)) & (data["Hugo_Symbol"].isin(gene_set)), :]
    data = data.loc[(data["Chromosome"].isin(step00.chromosome_list)), :]
    data[sample] = data["t_alt_count"] / data["t_depth"]
    data.set_index(keys=["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Hugo_Symbol", "Variant_Classification", "HGVSp_Short"], inplace=True, verify_integrity=True)
    data = data[[sample]]

    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data data CSV file", type=str)
    parser.add_argument("census", help="Census data data CSV file", type=str)
    parser.add_argument("driver", help="MutEnricher Fisher enrichment output", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.05)

    group_histology = parser.add_mutually_exclusive_group(required=True)
    group_histology.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_histology.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.census.endswith(".csv"):
        raise ValueError("Census must end with .CSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.p < 1):
        raise ValueError("P-values must be (0, 1)")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    patients = set(map(lambda x: step00.get_patient(x.split("/")[-1].split(".")[0]), args.input))

    if args.SQC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
        patients = set(filter(lambda x: step00.get_patient(x) in histology, patients))
    elif args.ADC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
        patients = set(filter(lambda x: step00.get_patient(x) in histology, patients))
    else:
        raise Exception("Something went wrong!!")
    print(sorted(patients))

    args.input = list(filter(lambda x: step00.get_patient(x.split("/")[-1].split(".")[0]) in patients, args.input))
    print(len(args.input))

    census_data = pandas.read_csv(args.census, index_col=0)
    gene_set = set(census_data.index)
    print(census_data)

    driver_data = pandas.read_csv(args.driver, sep="\t")
    driver_data = driver_data.loc[(driver_data["Gene"].isin(gene_set))]
    for column in step00.MutEnricher_pval_columns:
        driver_data = driver_data.loc[(driver_data[column] < args.p)]
    driver_data.sort_values(by="Fisher_pval", ascending=False, ignore_index=True, inplace=True)
    gene_set = set(driver_data["Gene"])
    print(driver_data)

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(read, args.input), axis="columns", join="outer", sort=True, verify_integrity=True).fillna(0).T
    print(input_data)

    tsne_data = pandas.DataFrame(data=sklearn.manifold.TSNE(perplexity=50, init="pca", verbose=1, random_state=42, method="exact", n_jobs=args.cpus).fit_transform(input_data), index=input_data.index, columns=["tSNE1", "tSNE2"], dtype=float)
    for column in list(tsne_data.columns):
        tsne_data[column] = sklearn.preprocessing.scale(tsne_data[column])
    tsne_data["Stage"] = list(map(step00.get_long_sample_type, list(tsne_data.index)))
    tsne_data["Patient"] = list(map(step00.get_patient, list(tsne_data.index)))
    order = list(filter(lambda x: x in set(tsne_data["Stage"]), step00.long_sample_type_list))
    print(tsne_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.scatterplot(data=tsne_data, x="tSNE1", y="tSNE2", hue="Stage", style="Stage", legend="brief", s=1000, hue_order=order)

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
