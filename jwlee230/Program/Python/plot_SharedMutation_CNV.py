"""
plot_SharedMutation_CNV.py: plot shared mutation with CNV segment
"""
import argparse
import multiprocessing
import tarfile
import typing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import tqdm
import step00

watching = ""
input_data = pandas.DataFrame()
CNV_data = pandas.DataFrame()
CGC_data = pandas.DataFrame()

lower_precancer_list: typing.List[str] = list()
higher_precancer_list: typing.List[str] = list()
lower_primary_list: typing.List[str] = list()
higher_primary_list: typing.List[str] = list()


def run(MSP: str, gene: str) -> typing.Optional[str]:
    chromosome = CGC_data.loc[gene, "Chromosome"]
    start = CGC_data.loc[gene, "Start position"]
    end = CGC_data.loc[gene, "End position"]

    selected_input_data = input_data.loc[(input_data["Hugo_Symbol"] == gene) & (input_data["Precancer"].isin(lower_precancer_list + higher_precancer_list))]
    if "SYN" not in MSP:
        selected_input_data = selected_input_data.loc[(selected_input_data[step00.nonsynonymous_column].isin(step00.nonsynonymous_mutations))]

    if selected_input_data.empty:
        return None

    selected_CNV_data = CNV_data.loc[(CNV_data["Sample"].isin(lower_precancer_list + lower_primary_list + higher_precancer_list + higher_primary_list)) & (CNV_data["chromosome"] == chromosome) & (start <= CNV_data["end"]) & (CNV_data["start"] <= end)]

    fig, axs = matplotlib.pyplot.subplots(figsize=(36, 18), nrows=2, sharex=True, sharey=True)

    for sample in lower_precancer_list:
        shared_location_list = list(selected_input_data.loc[(selected_input_data["Precancer"] == sample), "Start_Position"])

        for index, row in selected_CNV_data.loc[(selected_CNV_data["Sample"] == sample)].iterrows():
            axs[0].plot([row["start"], row["end"]], [row[watching], row[watching]], color="tab:pink", linewidth=10)
            axs[0].plot([location for location in shared_location_list], [row[watching] for location in shared_location_list], "r*", markersize=40, zorder=step00.small)

        for index, row in selected_CNV_data.loc[(selected_CNV_data["Sample"] == step00.get_paired_primary(sample))].iterrows():
            axs[0].plot([row["start"], row["end"]], [row[watching], row[watching]], color="tab:gray", linewidth=10)
            axs[0].plot([location for location in shared_location_list], [row[watching] for location in shared_location_list], "r*", markersize=40, zorder=step00.small)

    for sample in higher_precancer_list:
        shared_location_list = list(selected_input_data.loc[(selected_input_data["Precancer"] == sample), "Start_Position"])

        for index, row in selected_CNV_data.loc[(selected_CNV_data["Sample"] == sample)].iterrows():
            axs[1].plot([row["start"], row["end"]], [row[watching], row[watching]], color="tab:pink", linewidth=10)
            axs[1].plot([location for location in shared_location_list], [row[watching] for location in shared_location_list], "r*", markersize=40, zorder=step00.small)

        for index, row in selected_CNV_data.loc[(selected_CNV_data["Sample"] == step00.get_paired_primary(sample))].iterrows():
            axs[1].plot([row["start"], row["end"]], [row[watching], row[watching]], color="tab:gray", linewidth=10)
            axs[1].plot([location for location in shared_location_list], [row[watching] for location in shared_location_list], "r*", markersize=40, zorder=step00.small)

    axs[0].set_xticks([])
    axs[1].set_xticks([])
    axs[0].set_xlim((start, end))
    axs[1].set_xlim((start, end))
    axs[0].set_xlabel(f"{gene} ({chromosome}:{start}-{end})")
    axs[1].set_xlabel(f"{gene} ({chromosome}:{start}-{end})")
    axs[0].set_ylabel("CNV of MSP-low")
    axs[1].set_ylabel("CNV of MSP-high")
    matplotlib.pyplot.tight_layout()

    file_name = f"{MSP}-{gene}.pdf"
    fig.savefig(file_name)
    matplotlib.pyplot.close(fig)
    return file_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input shared mutation information TSV file", type=str)
    parser.add_argument("CNV", help="CNV segment.tsv file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("cgc", help="CGC gene CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--watching", help="Watching column name", type=str, required=True)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)
    parser.add_argument("--threshold", help="Threshold for gain/loss", type=float, default=0.2)
    parser.add_argument("--percentage", help="Percentage of patients to include", type=float, default=0.1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.CNV.endswith(".tsv"):
        raise ValueError("CNV must end with .TSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.cgc.endswith(".csv"):
        raise ValueError("CGC must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.threshold < 1):
        raise ValueError("Threshold must be (0, 1)")
    elif not (0.0 < args.percentage < 0.5):
        raise ValueError("Percentage must be (0.0, 0.5)!!")

    watching = args.watching

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(len(patients), patients)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data = input_data.loc[(input_data["Patient"].isin(patients))]
    print(input_data)

    CNV_data = pandas.read_csv(args.CNV, sep="\t", index_col=0)
    CNV_data = CNV_data.loc[(CNV_data["Patient"].isin(patients))]
    print(CNV_data)

    CGC_data = pandas.read_csv(args.cgc, index_col=0)
    CGC_data = CGC_data.loc[~(CGC_data["Genome Location"].str.contains(":-", regex=False))]
    CGC_data["Chromosome"] = list(map(lambda x: "chr" + x[:x.find(":")], CGC_data["Genome Location"]))
    CGC_data["Start position"] = list(map(lambda x: int(x[x.find(":") + 1:x.find("-")]), CGC_data["Genome Location"]))
    CGC_data["End position"] = list(map(lambda x: int(x[x.find("-") + 1:]), CGC_data["Genome Location"]))
    print(CGC_data)

    gene_list = sorted(set(input_data["Hugo_Symbol"]) & set(CGC_data.index))
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

        lower_primary_list = list(map(step00.get_paired_primary, lower_precancer_list))
        higher_primary_list = list(map(step00.get_paired_primary, higher_precancer_list))

        with multiprocessing.Pool(processes=args.cpus) as pool:
            figures += list(filter(None, list(pool.starmap(run, [(MSP, gene) for gene in gene_list]))))

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
