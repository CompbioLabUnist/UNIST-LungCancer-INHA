"""
make_input_revolver.py: make input file for revolver
"""
import argparse
import gtfparse
import multiprocessing
import pandas

gencode_data = pandas.DataFrame()
gene_symbol_set = set()


def change_position_gene_name(seqname: str, start: str, end: str) -> str:
    d = gencode_data.loc[(gencode_data["seqname"] == seqname) & (gencode_data["start"] <= start) & (end <= gencode_data["end"]), "gene_name"]
    if d.empty:
        return ""
    else:
        return d.to_numpy()[0]


def is_driver(symbol: str) -> bool:
    return symbol in gene_symbol_set


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input loci.TSV files", type=str)
    parser.add_argument("gencode", help="Gencode annotation GTF file", type=str)
    parser.add_argument("driver", help="Driver gene TSV file (not necessarily TSV)", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.gencode.endswith(".gtf"):
        raise ValueError("Gencode must end with .GTF!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    try:
        input_data = pandas.read_csv(args.input, sep="\t", usecols=["mutation_id", "sample_id", "cluster_id", "cellular_prevalence"])
    except pandas.errors.EmptyDataError:
        input_data = pandas.DataFrame(columns=["mutation_id", "sample_id", "cluster_id", "cellular_prevalence"])
    input_data["seqname"] = list(map(lambda x: x.split(":")[0], input_data["mutation_id"]))
    input_data["start"] = list(map(lambda x: int(x.split(":")[1]), input_data["mutation_id"]))
    input_data["end"] = list(map(lambda x: int(x.split(":")[2]), input_data["mutation_id"]))
    print(input_data)

    gencode_data = gtfparse.read_gtf(args.gencode, usecols=["seqname", "start", "end", "gene_name"]).drop_duplicates()
    print(gencode_data)

    driver_data = pandas.read_csv(args.driver, sep="\t")
    driver_data = driver_data.loc[(driver_data["Fisher_pval"] < 0.05)]
    print(list(driver_data.columns))
    print(driver_data)
    gene_symbol_set = set(driver_data["Gene"])

    with multiprocessing.Pool(args.cpus) as pool:
        input_data["gene_name"] = pool.starmap(change_position_gene_name, input_data[["seqname", "start", "end"]].to_numpy())
        input_data["is.driver"] = pool.map(is_driver, input_data["gene_name"])

    output_data = pandas.DataFrame()
    output_data["patientID"] = list(map(lambda x: x.split(".")[0], input_data["sample_id"]))
    output_data["variantID"] = input_data["mutation_id"]
    output_data["CCF"] = input_data["cellular_prevalence"]
    output_data["is.clonal"] = True
    output_data["is.driver"] = input_data["is.driver"]
    output_data["Misc"] = input_data["gene_name"]
    output_data["cluster"] = input_data["cluster_id"]

    print(output_data)
    output_data.to_csv(args.output, sep="\t", index=False)
