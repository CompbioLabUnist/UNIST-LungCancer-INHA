"""
select_DEGs.py: select up-regulated & down-regulated DEGs
"""
import argparse
import numpy
import pandas


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str)
    parser.add_argument("output", help="Output Basename", type=str)
    parser.add_argument("--padj", help="P-value threshold", type=float, default=0.05)
    parser.add_argument("--fold", help="Fold change threshold", type=float, default=2)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--up", help="Select UP-regulated genes", action="store_true", default=False)
    group.add_argument("--down", help="Select DOWN-regulated genes", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not (0 < args.padj < 1):
        raise ValueError("Padj must be [0, 1]!!")

    input_data = pandas.read_csv(args.input, sep="\t")

    if args.up:
        input_data = input_data.loc[(input_data["log2FoldChange"] >= numpy.log2(args.fold)) & (input_data["pvalue"] < args.padj) & (input_data["padj"] < args.padj)].sort_values(by="log2FoldChange", ascending=False)
    elif args.down:
        input_data = input_data.loc[(input_data["log2FoldChange"] <= -1 * numpy.log2(args.fold)) & (input_data["pvalue"] < args.padj) & (input_data["padj"] < args.padj)].sort_values(by="log2FoldChange", ascending=True)

    print(input_data)

    if input_data.empty:
        pandas.DataFrame(columns=["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"], index=["None"], data=[["", "", "", "", "", ""]]).to_csv(args.output + ".tsv", sep="\t", float_format="%.2e")
        pandas.DataFrame(columns=["gene", "log2FoldChange", "pvalue", "padj"], index=[0], data=[["None", "", "", ""]]).to_latex(args.output + ".tex", index=False, float_format="%.2e")
    else:
        input_data.to_csv(args.output + ".tsv", sep="\t", float_format="%.2e")
        input_data.iloc[:3, :].loc[:, ["log2FoldChange", "pvalue", "padj"]].to_latex(args.output + ".tex", float_format="%.2e")
