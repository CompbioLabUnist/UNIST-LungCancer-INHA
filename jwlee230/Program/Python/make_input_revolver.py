"""
make_input_revolver.py: make input file for revolver
"""
import argparse
import multiprocessing
import typing
import pandas
import step00

input_data = pandas.DataFrame()
driver_data = pandas.DataFrame()
gene_symbol_set = set()


def read_tsv(filename: str) -> pandas.DataFrame:
    try:
        input_data = pandas.read_csv(filename, sep="\t", usecols=["mutation_id", "sample_id", "cluster_id", "cellular_prevalence"])
    except pandas.errors.EmptyDataError:
        return pandas.DataFrame(columns=["mutation_id", "sample_id", "cluster_id", "cellular_prevalence", "seqname", "start", "end"])

    input_data["seqname"] = list(map(lambda x: x.split(":")[0], input_data["mutation_id"]))
    input_data["start"] = list(map(lambda x: int(x.split(":")[1]), input_data["mutation_id"]))
    input_data["end"] = list(map(lambda x: int(x.split(":")[2]), input_data["mutation_id"]))

    return input_data


def change_position_Gene(seqname: str, start: str, end: str) -> typing.List[str]:
    d = driver_data.loc[(driver_data["seqname"] == seqname) & (driver_data["start"] <= start) & (end <= driver_data["end"]), "Gene"]
    if d.empty:
        return [""]
    else:
        return list(set(d.to_numpy()))


def is_driver(symbol: str) -> bool:
    return symbol in gene_symbol_set


def get_CCF(mutation_id: str, Gene: str) -> str:
    data = input_data.loc[(input_data["mutation_id"] == mutation_id) & (input_data["Gene"] == Gene), "cellular_prevalence"]
    answer = ""
    for i, d in enumerate(data.to_numpy()):
        answer += "R{0:d}:{1:f};".format(i + 1, d)

    if answer.endswith(";"):
        answer = answer[:-1]

    return answer


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input loci.TSV files", type=str, nargs="+")
    parser.add_argument("driver", help="Driver gene TSV file (not necessarily TSV)", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.01)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.p < 1):
        raise ValueError("P must be (0, 1)")

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(read_tsv, args.input), axis="index", copy=False)
    print(input_data)

    driver_data = pandas.read_csv(args.driver, sep="\t")
    driver_data = driver_data.loc[(driver_data["Fisher_pval"] < args.p)]
    driver_data["seqname"] = list(map(lambda x: x.replace("-", ":").split(":")[0], driver_data["coordinates"]))
    driver_data["start"] = list(map(lambda x: int(x.replace("-", ":").split(":")[1]), driver_data["coordinates"]))
    driver_data["end"] = list(map(lambda x: int(x.replace("-", ":").split(":")[2]), driver_data["coordinates"]))
    print(list(driver_data.columns))
    print(driver_data)
    gene_symbol_set = set(driver_data["Gene"])

    with multiprocessing.Pool(args.cpus) as pool:
        input_data["Gene"] = pool.starmap(change_position_Gene, input_data[["seqname", "start", "end"]].to_numpy())
        input_data = input_data.explode("Gene", ignore_index=True)
        input_data["is.driver"] = pool.map(is_driver, input_data["Gene"])
    print(input_data)

    output_data = pandas.DataFrame()
    output_data["patientID"] = list(map(step00.get_patient, input_data["sample_id"]))
    output_data["variantID"] = list(map(lambda x: x[1] if x[0] == "" else x[0], input_data[["Gene", "mutation_id"]].to_numpy()))
    with multiprocessing.Pool(args.cpus) as pool:
        output_data["CCF"] = pool.starmap(get_CCF, input_data[["mutation_id", "Gene"]].to_numpy())
    output_data["is.clonal"] = "TRUE"
    output_data["is.driver"] = list(map(lambda x: "TRUE" if x else "FALSE", input_data["is.driver"]))
    output_data["Misc"] = input_data["mutation_id"]
    output_data["cluster"] = list(map(lambda x: "cluster{0:d}".format(x), input_data["cluster_id"]))
    output_data.drop_duplicates(inplace=True, ignore_index=True)

    print(output_data)
    output_data.to_csv(args.output, sep="\t", index=False)
