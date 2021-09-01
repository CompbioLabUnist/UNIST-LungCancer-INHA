"""
filter_fusioncatcher.py: Filtering FusionCatcher output
"""
import argparse
import pandas

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("normal", help="Normal sample TSV result", type=str)
    parser.add_argument("tumor", help="Tumor sample TSV result", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)

    args = parser.parse_args()

    if not args.normal.endswith(".tsv"):
        raise ValueError("Normal must end with .TSV!!")
    elif not args.tumor.endswith(".tsv"):
        raise ValueError("Tumor must end with .TSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")

    tumor_data = pandas.read_csv(args.tumor, sep="\t", index_col=["Gene_1_symbol(5end_fusion_partner)", "Gene_2_symbol(3end_fusion_partner)"])
    print(tumor_data)

    try:
        normal_data = pandas.read_csv(args.normal, sep="\t", index_col=["Gene_1_symbol(5end_fusion_partner)", "Gene_2_symbol(3end_fusion_partner)"])
        print(normal_data)
        filtered_data = tumor_data.drop(index=normal_data.index, errors="ignore")
    except pandas.errors.EmptyDataError:
        filtered_data = tumor_data
    print(filtered_data)

    filtered_data.to_csv(args.output, sep="\t")
