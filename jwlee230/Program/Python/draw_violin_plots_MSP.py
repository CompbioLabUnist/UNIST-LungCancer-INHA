"""
draw_violin_plots_MSP.py: draw violin plots upon DEG with MSP
"""
import argparse
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import statannotations.Annotator
import tqdm
import step00

input_data = pandas.DataFrame()


def run(MSP: str, gene: str) -> str:
    order = ["Lower", "Higher"]
    stage_order = list(filter(lambda x: all([(x in set(input_data.loc[(input_data[MSP] == compare), "Stage"])) for compare in order]), step00.long_sample_type_list))

    compare_list = list()
    for (c1, s1), (c2, s2) in [(("Lower", stage), ("Higher", stage)) for stage in stage_order] + [((compare, a), (compare, b)) for a, b in itertools.combinations(stage_order, r=2) for compare in order]:
        stat, p = scipy.stats.mannwhitneyu(input_data.loc[(input_data["Stage"] == s1) & (input_data[MSP] == c1), gene], input_data.loc[(input_data["Stage"] == s2) & (input_data[MSP] == c2), gene])
        if p < 0.05:
            compare_list.append(((c1, s1), (c2, s2)))

    try:
        stat, p = scipy.stats.kruskal(*[input_data.loc[(input_data["Stage"] == stage) & (input_data[MSP] == compare), gene] for compare in order for stage in stage_order])
    except ValueError:
        _, p = 0., 1.0

    if (p > 0.05) and (not compare_list):
        return ""

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.violinplot(data=input_data, x=MSP, order=order, y=gene, hue="Stage", hue_order=stage_order, palette=step00.stage_color_code, cut=1, linewidth=5, ax=ax)
    if compare_list:
        statannotations.Annotator.Annotator(ax, compare_list, data=input_data, x=MSP, order=order, y=gene, hue="Stage", hue_order=stage_order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.ylabel(f"{gene} expression")
    matplotlib.pyplot.title(f"{gene}: Kruskal-Wallis p={p:.3f}")
    matplotlib.pyplot.tight_layout()

    fig_name = f"{MSP}_{gene}.pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TPM TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("cgc", help="CGC CSV files", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_threshold = parser.add_mutually_exclusive_group(required=True)
    group_threshold.add_argument("--median", help="Median threshold", action="store_true", default=False)
    group_threshold.add_argument("--mean", help="Mean threshold", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .tsv!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.cgc.endswith(".csv"):
        raise ValueError("CGC must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

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

    cgc_data = pandas.read_csv(args.cgc, index_col="Gene Symbol")
    gene_set = set(cgc_data.index)
    print(cgc_data)

    input_data = pandas.read_csv(args.input, sep="\t", index_col="gene_name").T
    gene_set &= set(input_data.columns)
    input_data = input_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index))), sorted(gene_set)]
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    print(input_data)

    for MSP in tqdm.tqdm(step00.sharing_columns):
        if args.median:
            threshold = numpy.median(clinical_data[MSP])
        elif args.mean:
            threshold = numpy.mean(clinical_data[MSP])
        else:
            raise Exception("Something went wrong!!")
        print(MSP, threshold)

        input_data[MSP] = list(map(lambda x: "Lower" if (clinical_data.loc[step00.get_patient(x), MSP] < threshold) else "Higher", list(input_data.index)))

    with multiprocessing.Pool(args.cpus) as pool:
        figures = list(filter(None, pool.starmap(run, itertools.product(step00.sharing_columns, sorted(gene_set)))))

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(figures):
            tar.add(f, arcname=f)
