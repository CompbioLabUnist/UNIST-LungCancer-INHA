"""
compare_TIMER_sequenza.py: Compare TIMER & Sequenza results
"""
import argparse
import collections
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import tqdm
import step00

watching = "depth.ratio"
input_data = pandas.DataFrame()


def get_data(file_name: str) -> pandas.DataFrame:
    data = pandas.read_csv(file_name, sep="\t", usecols=["chromosome", "start.pos", "end.pos", watching]).dropna(axis="index")
    data["sample"] = file_name.split("/")[-2]
    return data


def run(cell: str) -> str:
    title, tool = cell.split("_")

    r, p = scipy.stats.pearsonr(input_data["CNV burden"], input_data[cell])

    try:
        g = seaborn.jointplot(data=input_data, x="CNV burden", y=cell, kind="reg", height=24, ratio=6)
    except numpy.core._exceptions._ArrayMemoryError:
        return ""
    g.fig.text(0.5, 0.75, "r={0:.3f}, p={1:.3f}".format(r, p), color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")
    g.plot_marginals(seaborn.histplot, kde=True, stat="probability", multiple="stack")
    g.set_axis_labels("CNV burden", f"Score by {tool}")
    g.fig.suptitle(title)
    g.fig.tight_layout()

    fig_name = cell.replace(" ", "").replace("(", "").replace(")", "").replace("/", "_") + ".pdf"
    g.savefig(fig_name)
    matplotlib.pyplot.close(g.fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="TIMER result CSV file", type=str)
    parser.add_argument("cnv", help="Sequenza result file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--ratio", help="Ratio threshold for CNV", type=float, default=0.2)

    args = parser.parse_args()

    if not args.input.endswith(".csv"):
        raise ValueError("Input must end with .csv!!")
    elif list(filter(lambda x: not x.endswith("sample_segments.txt"), args.cnv)):
        raise ValueError("CNV must end with sample_segments.txt!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .tar!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.ratio < 1):
        raise ValueError("ratio must be (0, 1)")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = pandas.read_csv(args.input, index_col=0).T
    cell_types = list(input_data.columns)
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    for stage in set(input_data["Stage"]):
        if len(input_data.loc[(input_data["Stage"] == stage)]) < 3:
            input_data = input_data.loc[~(input_data["Stage"] == stage)]
    print(input_data)

    with multiprocessing.Pool(args.cpus) as pool:
        sequenza_data = pandas.concat(objs=pool.map(get_data, args.cnv), axis="index", copy=False, ignore_index=True, verify_integrity=True)
    sequenza_data = sequenza_data.loc[(sequenza_data["sample"].isin(set(input_data.index))) & (sequenza_data["chromosome"].isin(step00.chromosome_list))]
    sequenza_data = sequenza_data.loc[(sequenza_data[watching] <= (1 - args.ratio)) | (sequenza_data[watching] >= (1 + args.ratio))]
    print(sequenza_data)

    counter: collections.Counter = collections.Counter(sequenza_data["sample"])
    input_data["CNV burden"] = list(map(lambda x: counter[x], list(input_data.index)))
    print(input_data)

    with multiprocessing.Pool(args.cpus) as pool:
        tar_files = list(filter(None, pool.map(run, cell_types)))

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(tar_files):
            tar.add(f, arcname=f)
