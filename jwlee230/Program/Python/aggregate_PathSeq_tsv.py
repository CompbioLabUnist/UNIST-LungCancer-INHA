"""
aggregate_PathSeq_tsv.py: aggregate PathSeq as tsv
"""
import argparse
import multiprocessing
import pandas
import tqdm
import step00

input_data = pandas.DataFrame()


def get_data(filename: str) -> pandas.DataFrame:
    data = pandas.read_csv(filename, sep="\t")
    data["ID"] = step00.get_id(filename)
    return data


def get_real_taxonomy(taxon: str) -> str:
    return taxon.split("|")[-1].replace("_", " ")


def query(sample: str, taxon: str) -> float:
    data = input_data.loc[(input_data["real_taxonomy"] == taxon) & (input_data["ID"] == sample), "score_normalized"]
    if data.empty:
        return 0.0
    else:
        return data.to_numpy()[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PathSeq results TSV file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--level", choices=step00.PathSeq_type_list, type=str, required=True)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    args.input.sort(key=step00.sorting)

    sample_list = list(map(step00.get_id, args.input))
    print(len(sample_list), sample_list)

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="index", copy=False, ignore_index=True, verify_integrity=True)
        input_data["real_taxonomy"] = pool.map(get_real_taxonomy, input_data["taxonomy"])
    input_data = input_data.loc[(input_data["kingdom"] == "Bacteria") & (input_data["type"] == args.level)]
    print(input_data)

    taxa_list = sorted(set(input_data["real_taxonomy"]))

    output_data = pandas.DataFrame(index=sample_list, columns=taxa_list, dtype=float)
    with multiprocessing.Pool(args.cpus) as pool:
        for index in tqdm.tqdm(sample_list):
            output_data.loc[index, :] = pool.starmap(query, [(index, taxon) for taxon in taxa_list])

    for index in tqdm.tqdm(sample_list):
        output_data.loc[index, :] = output_data.loc[index, :] / sum(output_data.loc[index, :]) * 100

    output_data["Subtype"] = list(map(step00.get_long_sample_type, list(output_data.index)))
    print(output_data)

    output_data.to_csv(args.output, sep="\t")
