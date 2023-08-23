"""
draw_volcano_plot.py: draw the volcano plot of DEG
"""
import argparse
from adjustText import adjust_text
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("DEG", help="DEG TSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--padj", help="P-value threshold", type=float, default=0.05)
    parser.add_argument("--fold", help="Fold change threshold", type=float, default=2)
    parser.add_argument("--annotation", help="Annotation limits", type=int, default=10)

    args = parser.parse_args()

    if not args.DEG.endswith(".tsv"):
        raise ValueError("DEG must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif args.annotation < 1:
        raise ValueError("Annotation must be positive!!")

    DEG_data = pandas.read_csv(args.DEG, sep="\t", index_col=0).dropna(axis="index", how="any")
    DEG_data["-log(Padj)"] = -1 * numpy.log10(DEG_data["padj"], dtype=float)
    DEG_data["importance"] = list(map(lambda x: abs(x[0] * x[1]), zip(DEG_data["log2FoldChange"], DEG_data["-log(Padj)"])))
    DEG_data.sort_values(by="importance", ascending=False, inplace=True)
    print(DEG_data)

    up_gene = DEG_data.loc[(DEG_data["log2FoldChange"] >= numpy.log2(args.fold)) & (DEG_data["padj"] < args.padj) & (DEG_data["pvalue"] < args.padj), ["log2FoldChange", "-log(Padj)"]]
    down_gene = DEG_data.loc[(DEG_data["log2FoldChange"] <= -1 * numpy.log2(args.fold)) & (DEG_data["padj"] < args.padj) & (DEG_data["pvalue"] < args.padj), ["log2FoldChange", "-log(Padj)"]]
    NS_gene = DEG_data.loc[((DEG_data["log2FoldChange"] < numpy.log2(args.fold)) & (DEG_data["log2FoldChange"] > -1 * numpy.log2(args.fold))) | (DEG_data["padj"] >= args.padj), ["log2FoldChange", "-log(Padj)"]]

    texts = list()

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    matplotlib.pyplot.scatter(NS_gene["log2FoldChange"], NS_gene["-log(Padj)"], color="tab:gray")
    matplotlib.pyplot.scatter(up_gene["log2FoldChange"], up_gene["-log(Padj)"], color="tab:red")
    matplotlib.pyplot.scatter(down_gene["log2FoldChange"], down_gene["-log(Padj)"], color="tab:blue")

    matplotlib.pyplot.axhline(y=-1 * numpy.log10(args.padj), linestyle="--", color="black")
    matplotlib.pyplot.text(x=0, y=-1 * numpy.log10(args.padj), s=f"Padj={args.padj:.2f}", horizontalalignment="center", verticalalignment="baseline", fontsize="xx-small")

    matplotlib.pyplot.axvline(x=numpy.log2(args.fold), linestyle="--", color="black")
    matplotlib.pyplot.text(x=numpy.log2(args.fold), y=-1 * numpy.log10(args.padj), s=f"log2(FC)={numpy.log2(args.fold):.1f}", rotation="vertical", horizontalalignment="left", verticalalignment="bottom", fontsize="xx-small")

    matplotlib.pyplot.axvline(x=-1 * numpy.log2(args.fold), linestyle="--", color="black")
    matplotlib.pyplot.text(x=-1 * numpy.log2(args.fold), y=-1 * numpy.log10(args.padj), s=f"log2(FC)={-1 * numpy.log2(args.fold):.1f}", rotation="vertical", horizontalalignment="right", verticalalignment="bottom", fontsize="xx-small")

    for index, d in tqdm.tqdm(up_gene.iloc[:args.annotation, :].iterrows()):
        texts.append(matplotlib.pyplot.text(s=index, x=d["log2FoldChange"], y=d["-log(Padj)"], color="tab:red", fontsize="large"))
    for index, d in tqdm.tqdm(down_gene.iloc[:args.annotation, :].iterrows()):
        texts.append(matplotlib.pyplot.text(s=index, x=d["log2FoldChange"], y=d["-log(Padj)"], color="tab:blue", fontsize="large"))

    adjust_text(texts, arrowprops={"arrowstyle": "-", "color": "k", "linewidth": 1, "alpha": 0.3}, ax=ax, lim=step00.small)

    limit = numpy.ceil(max(abs(DEG_data["log2FoldChange"])))

    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.xlabel("Normal-like← log2(FoldChange) →Primary-like")
    matplotlib.pyplot.ylabel("-log10(P.adj.)")
    matplotlib.pyplot.title("Up: {0:d} & Down: {1:d}".format(len(up_gene), len(down_gene)))
    matplotlib.pyplot.xlim(-1 * limit, limit)
    if matplotlib.pyplot.ylim()[1] < 2:
        matplotlib.pyplot.ylim(top=2)
    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
