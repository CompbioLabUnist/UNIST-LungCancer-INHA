"""
aggregate_RSEM_2.py: aggregate RSEM results for recurrence/non-recurrence
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
    data = pandas.read_csv(filename, sep="\t", usecols=["gene_id", "TPM"], index_col="gene_id")
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


def recur_or_non(recur, non, ID):
    recur = list(map(lambda x: x.split("/")[-1].split(".")[0], recur))
    non = list(map(lambda x: x.split("/")[-1].split(".")[0], non))
    if ID in recur:
        return "Recurrence"
    elif ID in non:
        return "Non-recurrence"
    else:
        raise Exception("Something went wrong!!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("gencode", help="Gencode annotation GTF file", type=str)
    parser.add_argument("trembl", help="Gencode TREMBL gz file", type=str)
    parser.add_argument("output", help="Output file basename", type=str)
    parser.add_argument("--recur", help="Recurrences RSEM genes.results file(s)", type=str, nargs="+", required=True)
    parser.add_argument("--non", help="Recurrences RSEM genes.results file(s)", type=str, nargs="+", required=True)
    parser.add_argument("--target", help="Give compared histologies", required=True)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if not args.gencode.endswith(".gtf"):
        raise ValueError("Gencode must end with .gtf!!")
    elif list(filter(lambda x: not x.endswith(".genes.results"), args.recur + args.non)):
        raise ValueError("Recur/Non-recur must end with .genes.results!!")
    elif args.target not in step00.long_sample_type_list:
        raise ValueError("TARGET is not valid!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    print(args.recur)
    args.recur = list(filter(lambda x: step00.get_long_sample_type(x.split("/")[-1].split(".")[0]) == args.target, args.recur))
    args.non = list(filter(lambda x: step00.get_long_sample_type(x.split("/")[-1].split(".")[0]) == args.target, args.non))
    input_files = sorted(args.non) + sorted(args.recur)
    print(input_files)

    trembl_data = pandas.read_csv(args.trembl, sep="\t")
    trembl_data = trembl_data.loc[(trembl_data["db_name"] == "Uniprot/SPTREMBL"), ["xref", "gene_stable_id"]]
    trembl_ID_set = set(trembl_data["xref"])
    print(trembl_data)

    gencode_data = gtfparse.read_gtf(args.gencode)[["gene_id", "gene_name"]].drop_duplicates(ignore_index=True)
    print(gencode_data)

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(read_RSEM, input_files), axis="columns", sort=True, verify_integrity=True, copy=False)
        input_data["ENSG"] = pool.map(trembl_to_ENSG, list(input_data.index))
        input_data = input_data.loc[(input_data["ENSG"] != "")]
        input_data["gene_name"] = pool.map(ENSG_to_names, list(input_data["ENSG"]))
        input_data = input_data.loc[(input_data["gene_name"] != "")]
    del input_data["ENSG"]
    print(input_data)
    print(list(input_data.columns))

    input_data.groupby(by="gene_name").sum().to_csv(args.output + ".tsv", sep="\t", index=True, header=True)
    pandas.DataFrame(data=[(ID, recur_or_non(args.recur, args.non, ID)) for ID in list(input_data.columns)[:-1]], columns=["ID", "condition"]).to_csv(args.output + ".coldata", sep="\t", index=False, header=True)
