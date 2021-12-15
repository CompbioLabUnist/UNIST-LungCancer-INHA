"""
draw_venn_plot_3.py: draw a venn diagram of DEG include every updown
"""
import argparse
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import venn
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("DEG", help="DEG TSV file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--annotation", help="Annotation for venn diagram", type=str, nargs="+", required=True)
    parser.add_argument("--padj", help="P-value threshold", type=float, default=0.05)
    parser.add_argument("--fold", help="Fold change threshold", type=float, default=2)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.DEG)):
        raise ValueError("DEG must end with .TSV!!")
    elif len(args.DEG) != len(args.annotation):
        raise ValueError("Annotation must be one-to-one upon DEG!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    input_data = dict()
    for annotation, input_file in zip(args.annotation, args.DEG):
        DEG_data = pandas.read_csv(input_file, sep="\t", header=0, names=["gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"], index_col="gene_id").dropna(axis="index", how="any")
        print(DEG_data)

        tmp_data = DEG_data.loc[(DEG_data["log2FoldChange"] >= numpy.log2(args.fold)) & (DEG_data["padj"] < args.padj) & (DEG_data["pvalue"] < args.padj), :]
        input_data[annotation + "_Up"] = set(tmp_data.index)
        print(tmp_data)

        tmp_data = DEG_data.loc[(DEG_data["log2FoldChange"] <= -1 * numpy.log2(args.fold)) & (DEG_data["padj"] < args.padj) & (DEG_data["pvalue"] < args.padj), :]
        input_data[annotation + "_Down"] = set(tmp_data.index)
        print(tmp_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))
    if len(args.DEG) == 3:
        venn.pseudovenn(input_data, ax=ax, fmt=step00.venn_format, fontsize=step00.matplotlib_parameters["legend.fontsize"], legend_loc="upper left")
    elif 0 < len(args.DEG) < 3:
        venn.venn(input_data, ax=ax, fmt=step00.venn_format, fontsize=step00.matplotlib_parameters["legend.fontsize"], legend_loc="upper left")
    else:
        raise Exception("Something went wrong!!")

    matplotlib.pyplot.tight_layout()
    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
