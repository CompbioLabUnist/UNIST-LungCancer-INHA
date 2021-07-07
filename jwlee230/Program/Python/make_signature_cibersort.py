"""
make_signature_cibersort.py: make the signature matrix file for CIBERSORTx
"""
import argparse
import collections
import pandas

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("annotation", help="Annotation TSV file (not necessarily TSV)", type=str)
    parser.add_argument("expression", help="Expression TSV file (not necessarily TSV)", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)

    args = parser.parse_args()

    if not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")

    annotation_data = pandas.read_csv(args.annotation, sep="\t", dtype=str, verbose=True, index_col="Index").dropna(axis="index", how="any")
    annotation_data = annotation_data.loc[(annotation_data["Sample_Origin"].isin(["nLung", "tLung"]))]
    print(annotation_data)

    for column in ["Sample_Origin", "Cell_type", "Cell_type.refined", "Cell_subtype"]:
        print(column, ":", collections.Counter(annotation_data[column]).most_common())

    expression_data = pandas.read_csv(args.expression, sep="\t", index_col="Index", usecols=["Index"] + list(annotation_data.index), verbose=True)
    expression_data = expression_data[annotation_data.index].T
    expression_data["Cell_subtype"] = annotation_data["Cell_subtype"]
    print(expression_data)

    grouped_data = expression_data.groupby(by="Cell_subtype").mean().T
    print(grouped_data)
    grouped_data.to_csv(args.output, sep="\t")
