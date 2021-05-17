"""
draw_volcano_plot.py: draw the volcano plot of DEG
"""
import argparse
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("DEG", help="DEG TSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--pvalue", help="P-value threshold", type=float, default=0.05)
    parser.add_argument("--fold", help="Fold change threshold", type=float, default=2)

    args = parser.parse_args()

    if not args.DEG.endswith(".tsv"):
        raise ValueError("DEG must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    DEG_data = pandas.read_csv(args.DEG, sep="\t", header=0, names=["gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"], index_col="gene_id")
    DEG_data["-logP"] = -1 * numpy.log10(DEG_data["pvalue"], dtype=float)
    print(DEG_data)

    up_gene = DEG_data.loc[(DEG_data["log2FoldChange"] >= numpy.log2(args.fold)) & (DEG_data["pvalue"] < args.pvalue), ["log2FoldChange", "-logP"]]
    down_gene = DEG_data.loc[(DEG_data["log2FoldChange"] <= -1 * numpy.log2(args.fold)) & (DEG_data["pvalue"] < args.pvalue), ["log2FoldChange", "-logP"]]
    NS_gene = DEG_data.loc[((DEG_data["log2FoldChange"] > -1 * numpy.log2(args.fold)) & (DEG_data["log2FoldChange"] < numpy.log2(args.fold))) | (DEG_data["pvalue"] >= args.pvalue), ["log2FoldChange", "-logP"]]

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
    matplotlib.pyplot.scatter(NS_gene["log2FoldChange"], NS_gene["-logP"], c="tab:gray")
    matplotlib.pyplot.scatter(up_gene["log2FoldChange"], up_gene["-logP"], c="tab:red")
    matplotlib.pyplot.scatter(down_gene["log2FoldChange"], down_gene["-logP"], c="tab:blue")

    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.xlabel("log2(Fold_Change)")
    matplotlib.pyplot.ylabel("-log10(P)")
    matplotlib.pyplot.title("Up: {0:d}, Down: {1:d}".format(len(up_gene), len(down_gene)))

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
