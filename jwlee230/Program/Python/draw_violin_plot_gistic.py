"""
draw_violin_plot_gistic.py: draw Violin plots with gistic result
"""
import argparse
import multiprocessing
import tarfile
import typing
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
stage_order = list()
state_palette = list()


def run(gene: str) -> typing.Tuple[str, float]:
    print("Running:", gene)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    gene_name, cytoband = gene.split(";")

    seaborn.violinplot(data=input_data, x="Stage", y=gene, order=stage_order, palette=stage_palette)
    statannotations.Annotator.Annotator(ax, list(zip(stage_order, stage_order[1:])), data=input_data, x="Stage", y=gene, order=stage_order).configure(test="Mann-Whitney", text_format="star", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.title("{0} in {1}".format(gene_name, cytoband))
    matplotlib.pyplot.ylabel("Segment Mean")
    matplotlib.pyplot.tight_layout()

    fig_name = "{1}_{0}.pdf".format(gene_name, cytoband)
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    final_p = 1.0
    for stage_a, stage_b in zip(stage_order, stage_order[1:]):
        stat, p = scipy.stats.mannwhitneyu(input_data.loc[(input_data["Stage"] == stage_a), gene], input_data.loc[(input_data["Stage"] == stage_b), gene])
        final_p = min(final_p, p)

    return (fig_name, final_p)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Gistic result TSV file (not necessarily TSV)", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.05)

    args = parser.parse_args()

    if not args.output.endswith(".tar"):
        raise ValueError("Output must end with .tar!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.p < 1):
        raise ValueError("P-value must be between 0 and 1!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = pandas.read_csv(args.input, sep="\t", index_col="Gene Symbol").T
    input_data.columns = gene_list = list(map(lambda x: "{0};{1}".format(x[0].split("|")[0], x[1]), zip(list(input_data.columns), input_data.loc["Cytoband", :])))
    input_data = input_data.iloc[2:, :]
    for column in tqdm.tqdm(gene_list):
        input_data[column] = numpy.power(2, input_data[column])
    input_data = input_data.astype(float)
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    print(input_data)

    stage_set = set(input_data["Stage"])
    stage_order = list(filter(lambda x: x in stage_set, step00.long_sample_type_list))
    stage_palette = list(map(lambda x: step00.stage_color_code[x], stage_order))

    with multiprocessing.Pool(args.cpus) as pool:
        files = list(filter(lambda x: x[1] < args.p, pool.map(run, gene_list)))

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(files):
            tar.add(f, arcname=f)
