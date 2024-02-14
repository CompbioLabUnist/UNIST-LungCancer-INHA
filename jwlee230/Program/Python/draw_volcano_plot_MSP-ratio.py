"""
draw_volcano_plot_MSP-ratio.py: draw the volcano plot of gene ratio & MSP
"""
import argparse
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
    for MSP in tqdm.tqdm(step00.sharing_columns):
        if f"{MSP}-importance" not in set(input_data.columns):
            continue

        input_data.sort_values(by=f"{MSP}-importance", ascending=False, inplace=True)

        up_gene = input_data.loc[(input_data[f"{MSP}-r"] > args.r) & (input_data[f"{MSP}-slope"] > args.slope)]
        down_gene = input_data.loc[(input_data[f"{MSP}-r"] < (-1 * args.r)) & (input_data[f"{MSP}-slope"] < (-1 * args.slope))]
        NS_gene = input_data.loc[list(filter(lambda x: (x not in set(up_gene)) and (x not in set(down_gene)), list(input_data.index)))]

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        matplotlib.pyplot.scatter(NS_gene[f"{MSP}-r"], NS_gene[f"{MSP}-log10(abs(slope))"], color="tab:gray")
        matplotlib.pyplot.scatter(up_gene[f"{MSP}-r"], up_gene[f"{MSP}-log10(abs(slope))"], color="tab:red")
        matplotlib.pyplot.scatter(down_gene[f"{MSP}-r"], down_gene[f"{MSP}-log10(abs(slope))"], color="tab:blue")

        matplotlib.pyplot.axhline(y=numpy.log10(args.slope), linestyle="--", color="black")
        matplotlib.pyplot.text(x=0, y=numpy.log10(args.slope), s=f"slope={args.slope:.1f}", horizontalalignment="center", verticalalignment="baseline", fontsize="xx-small", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.axvline(x=args.r, linestyle="--", color="black")
        matplotlib.pyplot.text(x=args.r, y=numpy.log10(args.slope), s=f"r={args.r:.1f}", rotation="vertical", horizontalalignment="left", verticalalignment="bottom", fontsize="xx-small", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.axvline(x=-1 * args.r, linestyle="--", color="black")
        matplotlib.pyplot.text(x=-1 * args.r, y=numpy.log10(args.slope), s=f"r={-1 * args.r:.1f}", rotation="vertical", horizontalalignment="right", verticalalignment="bottom", fontsize="xx-small", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.title(f"{len(down_gene)} neg. & {len(up_gene)} pos.")
        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.xlabel("Pearson correlation r")
        matplotlib.pyplot.ylabel("log10(|slope|)")
        matplotlib.pyplot.xlim(-1, 1)
        matplotlib.pyplot.tight_layout()

        figures.append(f"{MSP}-MT.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    input_data = input_data.loc[list(filter(lambda x: not x.startswith("MT-"), list(input_data.index)))]

    for MSP in tqdm.tqdm(step00.sharing_columns):
        if f"{MSP}-importance" not in set(input_data.columns):
            continue

        input_data.sort_values(by=f"{MSP}-importance", ascending=False, inplace=True)

        up_gene = input_data.loc[(input_data[f"{MSP}-r"] > args.r) & (input_data[f"{MSP}-slope"] > args.slope)]
        down_gene = input_data.loc[(input_data[f"{MSP}-r"] < (-1 * args.r)) & (input_data[f"{MSP}-slope"] < (-1 * args.slope))]
        NS_gene = input_data.loc[list(filter(lambda x: (x not in set(up_gene)) and (x not in set(down_gene)), list(input_data.index)))]

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        matplotlib.pyplot.scatter(NS_gene[f"{MSP}-r"], NS_gene[f"{MSP}-log10(abs(slope))"], color="tab:gray")
        matplotlib.pyplot.scatter(up_gene[f"{MSP}-r"], up_gene[f"{MSP}-log10(abs(slope))"], color="tab:red")
        matplotlib.pyplot.scatter(down_gene[f"{MSP}-r"], down_gene[f"{MSP}-log10(abs(slope))"], color="tab:blue")

        matplotlib.pyplot.axhline(y=numpy.log10(args.slope), linestyle="--", color="black")
        matplotlib.pyplot.text(x=0, y=numpy.log10(args.slope), s=f"slope={args.slope:.1f}", horizontalalignment="center", verticalalignment="baseline", fontsize="xx-small", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.axvline(x=args.r, linestyle="--", color="black")
        matplotlib.pyplot.text(x=args.r, y=numpy.log10(args.slope), s=f"r={args.r:.1f}", rotation="vertical", horizontalalignment="left", verticalalignment="bottom", fontsize="xx-small", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.axvline(x=-1 * args.r, linestyle="--", color="black")
        matplotlib.pyplot.text(x=-1 * args.r, y=numpy.log10(args.slope), s=f"r={-1 * args.r:.1f}", rotation="vertical", horizontalalignment="right", verticalalignment="bottom", fontsize="xx-small", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.title(f"{len(down_gene)} neg. & {len(up_gene)} pos.")
        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.xlabel("Pearson correlation r")
        matplotlib.pyplot.ylabel("log10(|slope|)")
        matplotlib.pyplot.xlim(-1, 1)
        matplotlib.pyplot.tight_layout()

        figures.append(f"{MSP}-noMT.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
