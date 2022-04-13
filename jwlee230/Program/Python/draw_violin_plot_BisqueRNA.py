"""
draw_violin_plot_BisqueRNA.py: draw violin plot from BisqueRNA result
"""
import argparse
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import scipy.stats
import seaborn
import statannotations.Annotator
import tqdm
import step00

input_data = pandas.DataFrame()


def run(cell: str) -> str:
    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    tmp = set(input_data["Stage"])
    order = list(filter(lambda x: x in tmp, step00.long_sample_type_list))
    palette = list(map(lambda x: step00.stage_color_code[x], order))

    try:
        stat, p = scipy.stats.kruskal(*[input_data.loc[(input_data["Stage"] == stage), cell] for stage in order])
    except ValueError:
        stat, p = 0.0, 1.0

    seaborn.violinplot(data=input_data, x="Stage", y=cell, order=order, palette=palette)
    statannotations.Annotator.Annotator(ax, list(itertools.combinations(order, 2)), data=input_data, x="Stage", y=cell, order=order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.title(f"{cell}: Kruskal-Wallis p={p:.3f}")
    matplotlib.pyplot.ylabel("Proportion")
    matplotlib.pyplot.tight_layout()

    fig_name = cell.replace(" ", "").replace("(", "").replace(")", "").replace("/", "_") + ".pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="BisqueRNA result TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .tsv!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .tar!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0).T
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    for stage in set(input_data["Stage"]):
        if len(input_data.loc[(input_data["Stage"] == stage)]) < 3:
            input_data = input_data.loc[~(input_data["Stage"] == stage)]
    print(input_data)

    with multiprocessing.Pool(args.cpus) as pool:
        tar_files = pool.map(run, list(input_data.columns)[:-1])

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(tar_files):
            tar.add(f, arcname=f)
