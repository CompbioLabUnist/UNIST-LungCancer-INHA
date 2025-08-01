"""
aggregate_RSEM_deconvolution.py: aggregate RSEM results for deconvolution tools
"""
import argparse
import gtfparse
import multiprocessing
import pandas
import step00

trembl_data = pandas.DataFrame()
gencode_data = pandas.DataFrame()
trembl_ID_set = set()


def read_RSEM(filename: str) -> pandas.DataFrame:
    data = pandas.read_csv(filename, sep="\t", usecols=["gene_id", "expected_count"], index_col="gene_id")
    data.columns = [filename.split("/")[-1].split(".")[0]]
    return data


def trembl_to_ENSG(gene_id: str) -> str:
    if gene_id in trembl_ID_set:
        return trembl_data.loc[(trembl_data["xref"] == gene_id), "gene_stable_id"].to_numpy()[0]
    else:
        return ""


def ENSG_to_names(gene_id: str) -> str:
    result = gencode_data.loc[(gencode_data["gene_id"].str.contains(gene_id)), "gene_name"]
    if result.empty:
        return ""
    else:
        return result.to_numpy()[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input genes.results (TSV) files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("gencode", help="Gencode annotation GTF file", type=str)
    parser.add_argument("trembl", help="Gencode TREMBL TSV.gz file", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_histology = parser.add_mutually_exclusive_group(required=True)
    group_histology.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_histology.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".genes.results"), args.input)):
        raise ValueError("INPUT must end with .genes.results!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.gencode.endswith(".gtf"):
        raise ValueError("Gencode must end with .gtf!!")
    elif not args.trembl.endswith(".tsv.gz"):
        raise ValueError("TREMBL must end with .tsv.gz!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    patients = list(map(lambda x: x.split("/")[-1].split(".")[0], args.input))
    if args.SQC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
        patients = list(filter(lambda x: step00.get_patient(x) in histology, patients))
    elif args.ADC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
        patients = list(filter(lambda x: step00.get_patient(x) in histology, patients))
    else:
        raise Exception("Something went wrong!!")
    print(patients)
    args.input = list(filter(lambda x: x.split("/")[-1].split(".")[0] in patients, args.input))

    trembl_data = pandas.read_csv(args.trembl, sep="\t")
    trembl_data = trembl_data.loc[(trembl_data["db_name"] == "Uniprot/SPTREMBL"), ["xref", "gene_stable_id"]]
    trembl_ID_set = set(trembl_data["xref"])
    print(trembl_data)

    gencode_data = gtfparse.read_gtf(args.gencode)[["gene_id", "gene_name"]].drop_duplicates(ignore_index=True)
    print(gencode_data)

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(read_RSEM, args.input), axis="columns", sort=True, verify_integrity=True, copy=False)
        input_data["ENSG"] = pool.map(trembl_to_ENSG, list(input_data.index))
        input_data = input_data.loc[(input_data["ENSG"] != "")]
        input_data["gene_name"] = pool.map(ENSG_to_names, list(input_data["ENSG"]))
        input_data = input_data.loc[(input_data["gene_name"] != "")]
    del input_data["ENSG"]
    print(input_data)
    print(sorted(input_data.columns))

    input_data.groupby("gene_name").sum().to_csv(args.output, sep="\t", index=True, header=True)
