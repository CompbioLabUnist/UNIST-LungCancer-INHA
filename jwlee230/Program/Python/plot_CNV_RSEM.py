"""
plot_CNV_RSEM.py: scatter plot with CNV and RSEM gene expression
"""
import argparse
import multiprocessing
import tarfile
import typing
import warnings
import matplotlib
import matplotlib.pyplot
import pandas
import scipy.stats
import seaborn
import tqdm
import step00

warnings.filterwarnings("ignore")

RNA_data = pandas.DataFrame()
stage_list: typing.List[str] = list()
palette: typing.Dict[str, str] = dict()


def scatter_all(gene: str) -> typing.Tuple[str, float, float]:
    fig_name = f"{gene}_All.pdf"
    r, p = scipy.stats.pearsonr(RNA_data[gene], RNA_data["Ploidy"])

    g = seaborn.jointplot(data=RNA_data, x=gene, y="Ploidy", kind="scatter", hue="Stage", hue_order=stage_list, palette=palette, height=24, ratio=5)
    g.fig.text(0.5, 0.75, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")
    g.savefig(fig_name)
    matplotlib.pyplot.close(g.fig)

    return fig_name, r, p


def scatter_stage(stage: str, gene: str) -> typing.Tuple[str, float, float]:
    fig_name = f"{gene}_{stage}.pdf"
    drawing_data = RNA_data.loc[(RNA_data["Stage"] == stage)]

    try:
        r, p = scipy.stats.pearsonr(drawing_data[gene], drawing_data["Ploidy"])
    except ValueError:
        return fig_name, 0.0, 1.0

    g = seaborn.jointplot(data=drawing_data, x=gene, y="Ploidy", kind="reg", color=palette[stage], height=24, ratio=5)
    g.fig.text(0.5, 0.75, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")
    g.savefig(fig_name)
    matplotlib.pyplot.close(g.fig)

    return fig_name, r, p


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("CNV", help="CNV ploidy.tsv file", type=str)
    parser.add_argument("RNA", help="RSEM RNA-seq gene expression TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data CNV file", type=str)
    parser.add_argument("figure", help="Output TAR file", type=str)
    parser.add_argument("table", help="Output TSV file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)
    parser.add_argument("--threshold", help="Threshold for p-value", type=float, default=0.05)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.CNV.endswith(".tsv"):
        raise ValueError("CNV must end with .TSV!!")
    elif not args.RNA.endswith(".tsv"):
        raise ValueError("RNA must end with .TSV!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.figure.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif not args.table.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.threshold < 1):
        raise ValueError("Threshold must be (0, 1)")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    CNV_data = pandas.read_csv(args.CNV, sep="\t", index_col=0)
    CNV_data["Patient"] = list(map(step00.get_patient, list(CNV_data.index)))
    CNV_data = CNV_data.loc[(CNV_data["Patient"].isin(patients))]
    print(CNV_data)

    RNA_data = pandas.read_csv(args.RNA, sep="\t", index_col=0).T
    gene_list = list(RNA_data.columns)[:100]
    RNA_data["Patient"] = list(map(step00.get_patient, list(RNA_data.index)))
    RNA_data = RNA_data.loc[(RNA_data["Patient"].isin(patients))]
    print(RNA_data)

    samples = sorted(set(CNV_data.index) & set(RNA_data.index), key=step00.sorting)
    CNV_data = CNV_data.loc[samples, :]
    RNA_data = RNA_data.loc[samples, :]
    RNA_data["Ploidy"] = CNV_data["Ploidy"]
    RNA_data["Stage"] = list(map(step00.get_long_sample_type, list(RNA_data.index)))
    print(RNA_data)

    stage_list = list(filter(lambda x: x in set(RNA_data["Stage"]), step00.long_sample_type_list))
    palette = dict([(stage, step00.stage_color_code[stage]) for stage in stage_list])
    print(stage_list)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    output_data = pandas.DataFrame(index=gene_list)
    with multiprocessing.Pool(args.cpus) as pool:
        output_data[["All-fig", "All-r", "All-p"]] = pool.map(scatter_all, gene_list)
        for stage in tqdm.tqdm(stage_list):
            output_data[[f"{stage}-fig", f"{stage}-r", f"{stage}-p"]] = pool.starmap(scatter_stage, [(stage, gene) for gene in gene_list])
    print(output_data)

    with tarfile.open(args.figure, "w") as tar:
        for figure in tqdm.tqdm(output_data.loc[(output_data["All-p"] < args.threshold), "All-fig"]):
            tar.add(figure, arcname=figure)
        for stage in stage_list:
            for figure in tqdm.tqdm(output_data.loc[(output_data[f"{stage}-p"] < args.threshold), f"{stage}-fig"]):
                tar.add(figure, arcname=figure)

    output_data.to_csv(args.table, sep="\t")
