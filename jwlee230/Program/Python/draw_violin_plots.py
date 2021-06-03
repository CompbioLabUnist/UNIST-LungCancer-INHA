"""
draw_violin_plots.py: draw violin plots upon DEG
"""
import argparse
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import scipy
import seaborn
import statannot
import step00

input_data = pandas.DataFrame()


def run(gene: str, ADC: bool = False, SQC: bool = False) -> str:
    print(gene)

    for stage_a, stage_b in itertools.combinations(set(input_data["Stage"]), 2):
        a = input_data.loc[(input_data["Stage"] == stage_a), gene]
        b = input_data.loc[(input_data["Stage"] == stage_b), gene]
        if a.empty or b.empty:
            return ""
        if scipy.stats.ttest_ind(a, b, equal_var=False)[1] >= 0.05:
            return ""

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    if ADC:
        seaborn.violinplot(data=input_data, x="Stage", y=gene, order=step00.ADC_stage_list, ax=ax)
        statannot.add_stat_annotation(ax, data=input_data, x="Stage", y=gene, order=step00.ADC_stage_list, test="t-test_ind", box_pairs=itertools.combinations(step00.ADC_stage_list, 2), text_format="star", loc="inside", verbose=0)
    elif SQC:
        seaborn.violinplot(data=input_data, x="Stage", y=gene, order=step00.SQC_stage_list, ax=ax)
        statannot.add_stat_annotation(ax, data=input_data, x="Stage", y=gene, order=step00.SQC_stage_list, test="t-test_ind", box_pairs=itertools.combinations(step00.SQC_stage_list, 2), text_format="star", loc="inside", verbose=0)
    else:
        raise Exception("Something went wrong!!")

    fig_name = gene + ".pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("count", help="Count TSV file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--pvalue", help="P-value threshold", type=float, default=0.05)
    parser.add_argument("--fold", help="Fold change threshold", type=float, default=2)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ADC", help="Draw ADC pathway", action="store_true", default=False)
    group.add_argument("--SQC", help="Draw SQC pathway", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.count)):
        raise ValueError("Count must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    for input_file in args.count:
        input_data = input_data.append(pandas.read_csv(input_file, sep="\t", index_col="gene_name").T)
    input_data.drop_duplicates(inplace=True)
    input_data = input_data[[c for c in list(input_data) if len(input_data[c].unique()) > 1]]
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    print(input_data)

    with multiprocessing.Pool(args.cpus) as pool:
        tar_files = pool.starmap(run, [(gene, args.ADC, args.SQC) for gene in list(input_data.columns)[:-1]])
        tar_files = list(filter(None, tar_files))

    with tarfile.open(args.output, "w") as tar:
        for f in tar_files:
            tar.add(f, arcname=f)
