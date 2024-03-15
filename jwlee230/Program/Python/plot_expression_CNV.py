"""
plot_expression_CNV.py: plot gene expression & CNV
"""
import argparse
import gtfparse
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
CNV_data = pandas.DataFrame()
DEG_data = pandas.DataFrame()
gene_data = pandas.DataFrame()

lower_precancer_list: typing.List[str] = list()
higher_precancer_list: typing.List[str] = list()
lower_primary_list: typing.List[str] = list()
higher_primary_list: typing.List[str] = list()


def query(sample: str, chromosome: str, start: int, end: int) -> float:
    length = end - start + 1
    a = list()
    weights = list()

    tmp_data = CNV_data.loc[(CNV_data["Sample"] == sample) & (CNV_data["chromosome"] == chromosome) & (CNV_data["start"] <= start) & (end <= CNV_data["end"]), :]
    for index, row in tmp_data.iterrows():
        a.append(row[args.watching])
        weights.append(1)

    tmp_data = CNV_data.loc[(CNV_data["Sample"] == sample) & (CNV_data["chromosome"] == chromosome) & (start <= CNV_data["start"]) & (CNV_data["end"] <= end), :]
    for index, row in tmp_data.iterrows():
        a.append(row[args.watching])
        weights.append((row["end"] - row["start"] + 1) / length)

    tmp_data = CNV_data.loc[(CNV_data["Sample"] == sample) & (CNV_data["chromosome"] == chromosome) & (CNV_data["start"] <= start) & (start <= CNV_data["end"]), :]
    for index, row in tmp_data.iterrows():
        a.append(row[args.watching])
        weights.append((row["end"] - start + 1) / length)

    tmp_data = CNV_data.loc[(CNV_data["Sample"] == sample) & (CNV_data["chromosome"] == chromosome) & (CNV_data["start"] <= end) & (end <= CNV_data["end"]), :]
    for index, row in tmp_data.iterrows():
        a.append(row[args.watching])
        weights.append((end - row["start"] + 1) / length)

    if a and weights:
        return numpy.average(a=a, weights=weights)
    else:
        return 1.0


def run_precancer(gene: str) -> str:
    selected_gene_data = gene_data[(gene_data["gene_id"] == gene)]
    chromosome = list(selected_gene_data["seqname"])[0]
    start = min(selected_gene_data["start"])
    end = max(selected_gene_data["end"])

    output_data = pandas.DataFrame(index=lower_precancer_list + higher_precancer_list)
    output_data["Expression"] = DEG_data.loc[output_data.index, gene]
    output_data["CNV"] = list(map(lambda x: query(x, chromosome, start, end), list(output_data.index)))

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    seaborn.regplot(data=output_data, x="CNV", y="Expression", scatter=True, fit_reg=True, color="tab:pink", ax=ax)

    matplotlib.pyplot.title(f"{gene} in Precancer")
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.tight_layout()

    fig_name = f"{gene}-Precancer.pdf"
    matplotlib.pyplot.savefig(fig_name)
    matplotlib.pyplot.close()

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("CNV", help="CNV segment.tsv file", type=str)
    parser.add_argument("DEG", help="DEG TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("gene", help="Gene GTF file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--watching", help="Watching column name", type=str, required=True)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)
    parser.add_argument("--threshold", help="Threshold for gain/loss", type=float, default=0.2)
    parser.add_argument("--percentage", help="Percentage of patients to include", type=float, default=0.25)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.CNV.endswith(".tsv"):
        raise ValueError("CNV must end with .TSV!!")
    elif not args.DEG.endswith(".tsv"):
        raise ValueError("DEG must end with .TSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.gene.endswith(".gtf"):
        raise ValueError("Gene must end with .GTF!!")
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
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(len(patients), patients)

    CNV_data = pandas.read_csv(args.CNV, sep="\t", index_col=0)
    CNV_data = CNV_data.loc[(CNV_data["Patient"].isin(patients))]
    print(CNV_data)

    DEG_data = pandas.read_csv(args.DEG, sep="\t", index_col=0).T
    DEG_data = DEG_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(DEG_data.index))), :]
    print(DEG_data)

    gene_data = gtfparse.read_gtf(args.gene)
    gene_data = gene_data.loc[(gene_data["feature"] == "exon")]
    print(gene_data)

    # gene_list = sorted(set(DEG_data.columns) & set(gene_data["gene_id"]))
    gene_list = ["ELOB", "H2AC4", "KRI1", "MZT2B", "PRDX2", "PTGES3", "TJP2", "ZNF148", "MAPK8IP3"]
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
            figures += list(filter(None, list(pool.map(run_precancer, gene_list))))

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
