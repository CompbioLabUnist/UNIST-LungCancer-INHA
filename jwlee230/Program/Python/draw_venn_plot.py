"""
draw_venn_plot.py: draw a venn diagram of DEG
"""
import argparse
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import upsetplot
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("DEG", help="DEG TSV file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--annotation", help="Annotation for venn diagram", type=str, nargs="+", required=True)
    parser.add_argument("--padj", help="P-value threshold", type=float, default=0.05)
    parser.add_argument("--fold", help="Fold change threshold", type=float, default=2)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--up", help="Draw ADC pathway", action="store_true", default=False)
    group.add_argument("--down", help="Draw SQC pathway", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.DEG)):
        raise ValueError("DEG must end with .TSV!!")
    elif len(args.DEG) != len(args.annotation):
        raise ValueError("Annotation must be one-to-one upon DEG!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    input_data = dict()
    for annotation, input_file in tqdm.tqdm(zip(args.annotation, args.DEG)):
        DEG_data = pandas.read_csv(input_file, sep="\t", index_col=0).dropna(axis="index", how="any")

        if args.up:
            DEG_data = DEG_data.loc[(DEG_data["log2FoldChange"] >= numpy.log2(args.fold)) & (DEG_data["padj"] < args.padj) & (DEG_data["pvalue"] < args.padj), :]
        elif args.down:
            DEG_data = DEG_data.loc[(DEG_data["log2FoldChange"] <= -1 * numpy.log2(args.fold)) & (DEG_data["padj"] < args.padj) & (DEG_data["pvalue"] < args.padj), :]
        else:
            raise Exception("Something went wrong!!")

        input_data[annotation] = set(DEG_data.index)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig = matplotlib.pyplot.figure(figsize=(10 * len(input_data) + 10, 24))

    try:
        upsetplot.plot(upsetplot.from_contents(input_data), fig=fig, show_counts="%d", show_percentages=True, element_size=None)
    except IndexError:
        pass

    fig.savefig(args.output, bbox_inches="tight")

    matplotlib.pyplot.close(fig)
