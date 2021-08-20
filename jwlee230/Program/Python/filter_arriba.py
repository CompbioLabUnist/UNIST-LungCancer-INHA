"""
filter_arriba.py: Filtering Arriba fusion.tsv results
"""
import argparse
import csv
import pandas

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("tumor", help="Tumor sample TSV result", type=str)
    parser.add_argument("normal", help="Normal sample TSV result", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)

    args = parser.parse_args()

    if not args.tumor.endswith(".tsv"):
        raise ValueError("Tumor must end with .TSV!!")
    elif not args.normal.endswith(".tsv"):
        raise ValueError("Normal must end with .TSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")

    wanted_columns = ["#gene1", "gene2", "breakpoint1", "breakpoint2"]

    tumor_data = pandas.read_csv(args.tumor, sep="\t")
    tumor_data = tumor_data.loc[(tumor_data["confidence"] == "high")]
    tumor_data["#gene1"] = list(map(lambda x: x.split(","), tumor_data["#gene1"]))
    tumor_data["gene2"] = list(map(lambda x: x.split(","), tumor_data["gene2"]))
    tumor_data = tumor_data.explode(column="#gene1", ignore_index=True).explode(column="gene2", ignore_index=True)
    tumor_data["#gene1"] = list(map(lambda x: x[:x.find("(")] if "(" in x else x, tumor_data["#gene1"]))
    tumor_data["gene2"] = list(map(lambda x: x[:x.find("(")] if "(" in x else x, tumor_data["gene2"]))
    print(tumor_data)

    try:
        normal_data = pandas.read_csv(args.normal, sep="\t")
        normal_data = normal_data.loc[(normal_data["confidence"] == "high")]
        normal_data["#gene1"] = list(map(lambda x: x.split(","), normal_data["#gene1"]))
        normal_data["gene2"] = list(map(lambda x: x.split(","), normal_data["gene2"]))
        normal_data = normal_data.explode(column="#gene1", ignore_index=True).explode(column="gene2", ignore_index=True)
        normal_data["#gene1"] = list(map(lambda x: x[:x.find("(")] if "(" in x else x, normal_data["#gene1"]))
        normal_data["gene2"] = list(map(lambda x: x[:x.find("(")] if "(" in x else x, normal_data["gene2"]))
    except FileNotFoundError:
        normal_data = pandas.DataFrame(columns=wanted_columns)
    print(normal_data)

    output_data = [("gene1", "gene2", "breakpoint1", "breakpoint2")] + sorted(set(tumor_data[wanted_columns].itertuples(index=False, name=None)) - set(normal_data[wanted_columns].itertuples(index=False, name=None)))
    print(output_data)

    with open(args.output, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(output_data)
