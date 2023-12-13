"""
draw_volcano_geneset_MSP.py: draw volcano with gene set on MSP
"""
import argparse
import itertools
import tarfile
import typing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input DEG-MSP TSV file", type=str)
    parser.add_argument("geneset", help="Gene set TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--percentage", help="Percentage of genes to draw", type=float, default=0.1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ADC", help="Draw ADC pathway", action="store_true", default=False)
    group.add_argument("--SQC", help="Draw SQC pathway", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.geneset.endswith(".tsv"):
        raise ValueError("Gene-set must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif not (0.0 < args.percentage < 0.5):
        raise ValueError("Percentage must be between 0.0 and 0.5!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)

    geneset_data = pandas.read_csv(args.geneset, sep="\t", index_col=0, header=0, names=["Index", "Contents"])
    print(geneset_data)

    notable_genes: typing.List[str] = sorted(list(filter(lambda x: x in set(input_data.index), geneset_data.loc["GENE_SYMBOLS", "Contents"].split(","))))
    print(len(notable_genes), notable_genes)

    input_data = input_data.loc[list(filter(lambda x: x in notable_genes, list(input_data.index)))]
    print(input_data)

    stages = step00.long_sample_type_list + ["Precancer", "All"]
    print(stages)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    figures = list()
    for stage, MSP in tqdm.tqdm(list(itertools.product(stages, step00.sharing_columns))):
        if f"{stage}-{MSP}-importance" not in set(input_data.columns):
            continue

        if stage == "All":
            color = "tab:blue"
        elif stage == "Precancer":
            color = "tab:pink"
        else:
            color = step00.stage_color_code[stage]

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        matplotlib.pyplot.scatter(input_data[f"{stage}-{MSP}-r"], input_data[f"{stage}-{MSP}-log10(abs(slope))"], color=color)

        median = numpy.median(input_data[f"{stage}-{MSP}-log10(abs(slope))"])
        matplotlib.pyplot.axhline(y=median, linestyle="--", color="black")
        matplotlib.pyplot.text(x=0, y=median, s=f"median(slope)={median:.3f}", horizontalalignment="center", verticalalignment="baseline", fontsize="xx-small", bbox={"alpha": 0.3, "color": "white"})

        percentage = numpy.quantile(input_data[f"{stage}-{MSP}-log10(abs(slope))"], 1 - args.percentage)
        matplotlib.pyplot.axhline(y=percentage, linestyle="--", color="black")
        matplotlib.pyplot.text(x=0, y=percentage, s=f"upper(slope)={percentage:.3f}", horizontalalignment="center", verticalalignment="baseline", fontsize="xx-small", bbox={"alpha": 0.3, "color": "white"})

        minimum = min(input_data[f"{stage}-{MSP}-log10(abs(slope))"])
        positive = list(filter(lambda x: x > 0, input_data[f"{stage}-{MSP}-r"]))

        median = numpy.median(positive) if positive else 0.0
        matplotlib.pyplot.axvline(x=median, linestyle="--", color="black")
        matplotlib.pyplot.text(x=median, y=minimum, s=f"median(r)={median:.3f}", rotation="vertical", horizontalalignment="left", verticalalignment="bottom", fontsize="xx-small", bbox={"alpha": 0.3, "color": "white"})

        percentage = numpy.quantile(positive, 1 - args.percentage) if positive else 0.0
        matplotlib.pyplot.axvline(x=percentage, linestyle="--", color="black")
        matplotlib.pyplot.text(x=percentage, y=minimum, s=f"upper(r)={percentage:.3f}", rotation="vertical", horizontalalignment="left", verticalalignment="bottom", fontsize="xx-small", bbox={"alpha": 0.3, "color": "white"})

        negative = list(filter(lambda x: x < 0, input_data[f"{stage}-{MSP}-r"]))

        median = numpy.median(negative) if negative else 0.0
        matplotlib.pyplot.axvline(x=median, linestyle="--", color="black")
        matplotlib.pyplot.text(x=median, y=minimum, s=f"median(r)={median:.3f}", rotation="vertical", horizontalalignment="left", verticalalignment="bottom", fontsize="xx-small", bbox={"alpha": 0.3, "color": "white"})

        percentage = numpy.quantile(negative, args.percentage) if negative else 0.0
        matplotlib.pyplot.axvline(x=percentage, linestyle="--", color="black")
        matplotlib.pyplot.text(x=percentage, y=minimum, s=f"lower(r)={percentage:.3f}", rotation="vertical", horizontalalignment="left", verticalalignment="bottom", fontsize="xx-small", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.xlabel("Pearson correlation r")
        matplotlib.pyplot.ylabel("log10(|slope|)")
        matplotlib.pyplot.xlim(-1, 1)
        matplotlib.pyplot.title(f"{len(input_data)} genes")
        matplotlib.pyplot.tight_layout()

        figures.append(f"{stage}-{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
