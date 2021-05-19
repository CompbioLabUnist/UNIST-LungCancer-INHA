"""
select_gene_MesKit.py: select significantly different genes for MesKit input
"""
import argparse
import pandas


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input Gene file", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)

    args = parser.parse_args()

    if not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")

    data = pandas.read_csv(args.input, sep="\t")
    data = data.loc[(data["P"] < 0.05), "GeneSymbol"]
    print(data)
    data.to_csv(args.output, sep="\t", header=False, index=False)
