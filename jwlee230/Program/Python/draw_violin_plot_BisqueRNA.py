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
import seaborn
import statannot
import step00

input_data = pandas.DataFrame()


def run(cell: str) -> str:
    print("Running:", cell)
    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
    tmp = set(input_data["Stage"])

    order = list(filter(lambda x: x in tmp, step00.long_sample_type_list))

    seaborn.violinplot(data=input_data, x="Stage", y=cell, order=order)
    statannot.add_stat_annotation(ax, data=input_data, x="Stage", y=cell, order=order, test="Mann-Whitney", box_pairs=itertools.combinations(order, 2), text_format="star", loc="inside", verbose=0)

    matplotlib.pyplot.title(cell)
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

    input_data = pandas.read_csv(args.input, sep="\t")
    input_data = input_data.set_index(list(input_data.columns)[0]).T
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    for stage in set(input_data["Stage"]):
        if len(input_data.loc[(input_data["Stage"] == stage)]) < 3:
            input_data = input_data.loc[~(input_data["Stage"] == stage)]
    print(input_data)

    with multiprocessing.Pool(args.cpus) as pool:
        tar_files = sorted(pool.map(run, list(input_data.columns)[:-1]))

    with tarfile.open(args.output, "w") as tar:
        for f in tar_files:
            tar.add(f, arcname=f)
