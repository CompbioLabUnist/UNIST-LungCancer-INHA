"""
make_input_revolver.py: make input file for revolver
"""
import argparse
import multiprocessing
import typing
import numpy
import pandas
import step00

input_data = pandas.DataFrame()
driver_data = pandas.DataFrame()
gene_symbol_set = set()
CCF_dict: typing.Any = dict()


def read_tsv(filename: str) -> pandas.DataFrame:
    try:
        input_data = pandas.read_csv(filename, sep="\t", usecols=["mutation_id", "sample_id", "cluster_id", "cellular_prevalence"])
    except pandas.errors.EmptyDataError:
        return pandas.DataFrame(columns=["mutation_id", "sample_id", "cluster_id", "cellular_prevalence", "seqname", "start", "end", "compared"])

    input_data["seqname"] = list(map(lambda x: x.split(":")[0], input_data["mutation_id"]))
    input_data["start"] = list(map(lambda x: int(x.split(":")[1]), input_data["mutation_id"]))
    input_data["end"] = list(map(lambda x: int(x.split(":")[2]), input_data["mutation_id"]))
    input_data["compared"] = filename.split("/")[-2]

    return input_data


def change_position_Gene(seqname: str, start: str, end: str) -> typing.List[str]:
    d = driver_data.loc[(driver_data["seqname"] == seqname) & (driver_data["start"] <= start) & (end <= driver_data["end"]), "Gene"]
    if d.empty:
        return [""]
    else:
        return list(set(d.to_numpy()))


def is_driver(symbol: str) -> bool:
    return symbol in gene_symbol_set


def get_CCF(sample_id: str, Gene: str) -> str:
    if step00.get_long_sample_type(sample_id) == "Primary":
        return ""
    assert CCF_dict[sample_id][Gene][0] and CCF_dict[sample_id][Gene][1]
    return "R1:{0:f};R2:{1:f}".format(numpy.mean(CCF_dict[sample_id][Gene][0]), numpy.mean(CCF_dict[sample_id][Gene][1]))


def get_misc(sample_id: str, Gene: str) -> str:
    return "/".join(list(input_data.loc[(input_data["sample_id"] == sample_id) & (input_data["Gene"] == Gene), "mutation_id"].to_numpy()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input loci.TSV files", type=str, nargs="+")
    parser.add_argument("driver", help="Driver gene TSV file (not necessarily TSV)", type=str)
    parser.add_argument("census", help="Cancer gene census CSV file", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.01)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.census.endswith(".csv"):
        raise ValueError("Census must end with .CSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.p < 1):
        raise ValueError("P must be (0, 1)")

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(read_tsv, args.input), axis="index", copy=False)
        input_data["sample_type"] = pool.map(step00.get_long_sample_type, input_data["sample_id"])
    input_data["sample_id"] = list(map(lambda x: x.split(".")[0], input_data["sample_id"]))
    print(list(input_data.columns))
    print(input_data)

    driver_data = pandas.read_csv(args.driver, sep="\t")
    driver_data = driver_data.loc[(driver_data["Fisher_pval"] < args.p)]
    driver_data["seqname"] = list(map(lambda x: x.replace("-", ":").split(":")[0], driver_data["coordinates"]))
    driver_data["start"] = list(map(lambda x: int(x.replace("-", ":").split(":")[1]), driver_data["coordinates"]))
    driver_data["end"] = list(map(lambda x: int(x.replace("-", ":").split(":")[2]), driver_data["coordinates"]))
    print(list(driver_data.columns))
    print(driver_data)

    gene_symbol_set = set(driver_data["Gene"])
    print("Gene set:", len(gene_symbol_set))

    census_data = pandas.read_csv(args.census)
    gene_symbol_set &= set(census_data["Gene Symbol"])
    print("Gene set:", len(gene_symbol_set))

    with multiprocessing.Pool(args.cpus) as pool:
        input_data["Gene"] = pool.starmap(change_position_Gene, input_data[["seqname", "start", "end"]].to_numpy())
        input_data = input_data.explode("Gene", ignore_index=True)
        input_data["is.driver"] = pool.map(is_driver, input_data["Gene"])
    print(input_data)

    for index, row in input_data.iterrows():
        sample_id, Gene, sample_type = row["sample_id"], row["Gene"], row["sample_type"]

        if sample_type == "Primary":
            continue

        if sample_id not in CCF_dict:
            CCF_dict[sample_id] = dict()

        if Gene not in CCF_dict[sample_id]:
            CCF_dict[sample_id][Gene] = [[], []]

        CCF_dict[sample_id][Gene][0].append(row["cellular_prevalence"])

        tmp_data = list(input_data.loc[(input_data["sample_id"] == step00.get_paired_primary(sample_id)) & (input_data["compared"] == sample_id) & (input_data["Gene"] == Gene) & (input_data["mutation_id"] == row["mutation_id"]), "cellular_prevalence"].to_numpy())
        assert len(tmp_data) == 1
        CCF_dict[sample_id][Gene][1] += tmp_data

    output_data = pandas.DataFrame()
    output_data["patientID"] = input_data["sample_id"]
    output_data["variantID"] = list(map(lambda x: x[1] if x[0] == "" else x[0], input_data[["Gene", "mutation_id"]].to_numpy()))
    with multiprocessing.Pool(args.cpus) as pool:
        output_data["CCF"] = pool.starmap(get_CCF, input_data[["sample_id", "Gene"]].to_numpy())
    output_data["is.clonal"] = "TRUE"
    output_data["is.driver"] = list(map(lambda x: "TRUE" if x else "FALSE", input_data["is.driver"]))
    with multiprocessing.Pool(args.cpus) as pool:
        output_data["Misc"] = pool.starmap(get_misc, input_data[["sample_id", "Gene"]].to_numpy())
    output_data["cluster"] = list(map(lambda x: "cluster{0:d}".format(x), input_data["cluster_id"]))
    output_data = output_data.loc[(input_data["sample_type"] != "Primary")]
    output_data.drop_duplicates(subset=["patientID", "variantID", "CCF"], inplace=True, ignore_index=True)

    print(output_data)
    output_data.to_csv(args.output, sep="\t", index=False)
