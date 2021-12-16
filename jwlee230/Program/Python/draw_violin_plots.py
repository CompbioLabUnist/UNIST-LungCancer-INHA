"""
draw_violin_plots.py: draw violin plots upon DEG
"""
import argparse
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy
import seaborn
import statannotations.Annotator
import tqdm
import step00

input_data = pandas.DataFrame()


def read_coldata(filename: str):
    return pandas.read_csv(filename, sep="\t")


def run(gene: str) -> str:
    print(gene)

    for stage in set(input_data["Stage"]):
        if numpy.var(input_data.loc[(input_data["Stage"] == stage), gene]) == 0:
            return ""

    for stage_a, stage_b in itertools.combinations(set(input_data["Stage"]), 2):
        if scipy.stats.mannwhitneyu(input_data.loc[(input_data["Stage"] == stage_a), gene], input_data.loc[(input_data["Stage"] == stage_b), gene])[1] < 0.001:
            break
    else:
        return ""

    stage_order = list(filter(lambda x: x in set(input_data["Stage"]), step00.long_sample_type_list))

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.violinplot(data=input_data, x="Stage", y=gene, order=stage_order, ax=ax)
    statannotations.Annotator.Annotator(ax, list(itertools.combinations(stage_order, 2)), data=input_data, x="Stage", y=gene, order=stage_order).configure(test="Mann-Whitney", text_format="star", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.tight_layout()

    fig_name = gene + ".pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TPM TSV file", type=str)
    parser.add_argument("coldata", help="Coldata file for selecting samples", type=str, nargs="+")
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .tsv!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = pandas.read_csv(args.input, sep="\t", index_col="gene_name")
    print(input_data)

    with multiprocessing.Pool(args.cpus) as pool:
        patient_data = pool.map(read_coldata, args.coldata)
    print(patient_data)

    with multiprocessing.Pool(args.cpus) as pool:
        tar_files = []

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(tar_files):
            tar.add(f, arcname=f)
