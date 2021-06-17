"""
draw_venn_plot_2.py: draw a venn diagram of DEG with pairwise
"""
import argparse
import matplotlib
import matplotlib.pyplot
import pandas
import venn
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("DEG", help="DEG TSV file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.05)
    parser.add_argument("--fold", help="Fold change threshold", type=float, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--up", help="Draw ADC pathway", action="store_true", default=False)
    group.add_argument("--down", help="Draw SQC pathway", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.DEG)):
        raise ValueError("DEG must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    input_data = dict()
    for input_file in args.DEG:
        DEG_data = pandas.read_csv(input_file, sep="\t", header=0, names=["gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"], index_col="gene_id").dropna(axis="index", how="any")
        print(DEG_data)

        if args.up:
            DEG_data = DEG_data.loc[(DEG_data["log2FoldChange"] > args.fold) & (DEG_data["pvalue"] < args.p) & (DEG_data["padj"] < args.p), :]
        elif args.down:
            DEG_data = DEG_data.loc[(DEG_data["log2FoldChange"] < (-1 * args.fold)) & (DEG_data["pvalue"] < args.p) & (DEG_data["padj"] < args.p), :]
        else:
            raise Exception("Something went wrong!!")

        print(DEG_data)
        input_data[input_file.split("/")[-1].split(".")[2]] = set(DEG_data.index)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
    venn.venn(input_data, ax=ax, fmt="{size:d}", fontsize=step00.matplotlib_parameters["font.size"])

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
