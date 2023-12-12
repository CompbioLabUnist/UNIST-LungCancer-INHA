"""
draw_volcano_plot_MSP.py: draw the volcano plot of gene & MSP
"""
import argparse
import itertools
import tarfile
from adjustText import adjust_text
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input DEG-MSP TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--r", help="r-value threshold", type=float, default=0.3)
    parser.add_argument("--slope", help="Slope threshold", type=float, default=5)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif not (0 < args.r < 1):
        raise ValueError("r-value must be in (0, 1)!!")
    elif args.slope <= 0:
        raise ValueError("Slope must be positive!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)

    stages = step00.long_sample_type_list + ["Precancer", "All"]

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    figures = list()
    for stage, MSP in tqdm.tqdm(list(itertools.product(stages, step00.sharing_columns))):
        if f"{stage}-{MSP}-importance" not in set(input_data.columns):
            continue

        input_data.sort_values(by=f"{stage}-{MSP}-importance", ascending=False, inplace=True)

        up_gene = input_data.loc[(input_data[f"{stage}-{MSP}-r"] > args.r) & (input_data[f"{stage}-{MSP}-slope"] > args.slope)]
        down_gene = input_data.loc[(input_data[f"{stage}-{MSP}-r"] < (-1 * args.r)) & (input_data[f"{stage}-{MSP}-slope"] > args.slope)]
        NS_gene = input_data.loc[(((-1 * args.r) < input_data[f"{stage}-{MSP}-r"]) & (input_data[f"{stage}-{MSP}-r"] < args.r)) | (input_data[f"{stage}-{MSP}-slope"] <= args.slope)]

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        matplotlib.pyplot.scatter(NS_gene[f"{stage}-{MSP}-r"], NS_gene[f"{stage}-{MSP}-log10(abs(slope))"], color="tab:gray")
        matplotlib.pyplot.scatter(up_gene[f"{stage}-{MSP}-r"], up_gene[f"{stage}-{MSP}-log10(abs(slope))"], color="tab:red")
        matplotlib.pyplot.scatter(down_gene[f"{stage}-{MSP}-r"], down_gene[f"{stage}-{MSP}-log10(abs(slope))"], color="tab:blue")

        matplotlib.pyplot.axhline(y=numpy.log10(args.slope), linestyle="--", color="black")
        matplotlib.pyplot.text(x=0, y=numpy.log10(args.slope), s=f"slope={args.slope:.1f}", horizontalalignment="center", verticalalignment="baseline", fontsize="xx-small", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.axvline(x=args.r, linestyle="--", color="black")
        matplotlib.pyplot.text(x=args.r, y=numpy.log10(args.slope), s=f"r={args.r:.1f}", rotation="vertical", horizontalalignment="left", verticalalignment="bottom", fontsize="xx-small", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.axvline(x=-1 * args.r, linestyle="--", color="black")
        matplotlib.pyplot.text(x=-1 * args.r, y=numpy.log10(args.slope), s=f"r={-1 * args.r:.1f}", rotation="vertical", horizontalalignment="left", verticalalignment="bottom", fontsize="xx-small", bbox={"alpha": 0.3, "color": "white"})

        texts = list()

        for index, d in up_gene.iloc[:5, :].iterrows():
            texts.append(matplotlib.pyplot.text(s=index, x=d[f"{stage}-{MSP}-r"], y=d[f"{stage}-{MSP}-log10(abs(slope))"], color="tab:red", horizontalalignment="right", fontsize="medium", bbox={"alpha": 0.3, "color": "white"}))
        for index, d in down_gene.iloc[:5, :].iterrows():
            texts.append(matplotlib.pyplot.text(s=index, x=d[f"{stage}-{MSP}-r"], y=d[f"{stage}-{MSP}-log10(abs(slope))"], color="tab:blue", horizontalalignment="left", fontsize="medium", bbox={"alpha": 0.3, "color": "white"}))

        adjust_text(texts, arrowprops={"arrowstyle": "-", "color": "k", "linewidth": 1, "alpha": 0.3}, ax=ax, lim=step00.small)

        matplotlib.pyplot.title(f"{stage}: {len(down_gene)} neg. & {len(up_gene)} pos.")
        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.xlabel("Pearson correlation r")
        matplotlib.pyplot.ylabel("log10(|slope|)")
        matplotlib.pyplot.xlim(-1, 1)
        matplotlib.pyplot.tight_layout()

        figures.append(f"{stage}-{MSP}-MT.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    input_data = input_data.loc[list(filter(lambda x: not x.startswith("MT-"), list(input_data.index)))]

    for stage, MSP in tqdm.tqdm(list(itertools.product(stages, step00.sharing_columns))):
        if f"{stage}-{MSP}-importance" not in set(input_data.columns):
            continue

        input_data.sort_values(by=f"{stage}-{MSP}-importance", ascending=False, inplace=True)

        up_gene = input_data.loc[(input_data[f"{stage}-{MSP}-r"] > args.r) & (input_data[f"{stage}-{MSP}-slope"] > args.slope)]
        down_gene = input_data.loc[(input_data[f"{stage}-{MSP}-r"] < (-1 * args.r)) & (input_data[f"{stage}-{MSP}-slope"] > args.slope)]
        NS_gene = input_data.loc[(((-1 * args.r) < input_data[f"{stage}-{MSP}-r"]) & (input_data[f"{stage}-{MSP}-r"] < args.r)) | (input_data[f"{stage}-{MSP}-slope"] <= args.slope)]

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        matplotlib.pyplot.scatter(NS_gene[f"{stage}-{MSP}-r"], NS_gene[f"{stage}-{MSP}-log10(abs(slope))"], color="tab:gray")
        matplotlib.pyplot.scatter(up_gene[f"{stage}-{MSP}-r"], up_gene[f"{stage}-{MSP}-log10(abs(slope))"], color="tab:red")
        matplotlib.pyplot.scatter(down_gene[f"{stage}-{MSP}-r"], down_gene[f"{stage}-{MSP}-log10(abs(slope))"], color="tab:blue")

        matplotlib.pyplot.axhline(y=numpy.log10(args.slope), linestyle="--", color="black")
        matplotlib.pyplot.text(x=0, y=numpy.log10(args.slope), s=f"slope={args.slope:.1f}", horizontalalignment="center", verticalalignment="baseline", fontsize="xx-small", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.axvline(x=args.r, linestyle="--", color="black")
        matplotlib.pyplot.text(x=args.r, y=numpy.log10(args.slope), s=f"r={args.r:.1f}", rotation="vertical", horizontalalignment="left", verticalalignment="bottom", fontsize="xx-small", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.axvline(x=-1 * args.r, linestyle="--", color="black")
        matplotlib.pyplot.text(x=-1 * args.r, y=numpy.log10(args.slope), s=f"r={-1 * args.r:.1f}", rotation="vertical", horizontalalignment="right", verticalalignment="bottom", fontsize="xx-small", bbox={"alpha": 0.3, "color": "white"})

        texts = list()

        for index, d in up_gene.iloc[:5, :].iterrows():
            texts.append(matplotlib.pyplot.text(s=index, x=d[f"{stage}-{MSP}-r"], y=d[f"{stage}-{MSP}-log10(abs(slope))"], color="tab:red", horizontalalignment="right", fontsize="medium", bbox={"alpha": 0.3, "color": "white"}))
        for index, d in down_gene.iloc[:5, :].iterrows():
            texts.append(matplotlib.pyplot.text(s=index, x=d[f"{stage}-{MSP}-r"], y=d[f"{stage}-{MSP}-log10(abs(slope))"], color="tab:blue", horizontalalignment="left", fontsize="medium", bbox={"alpha": 0.3, "color": "white"}))

        adjust_text(texts, arrowprops={"arrowstyle": "-", "color": "k", "linewidth": 1, "alpha": 0.3}, ax=ax, lim=step00.small)

        matplotlib.pyplot.title(f"{stage}: {len(down_gene)} neg. & {len(up_gene)} pos.")
        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.xlabel("Pearson correlation r")
        matplotlib.pyplot.ylabel("log10(|slope|)")
        matplotlib.pyplot.xlim(-1, 1)
        matplotlib.pyplot.tight_layout()

        figures.append(f"{stage}-{MSP}-noMT.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
