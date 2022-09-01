"""
merge_sequenza.py: merge Sequenza alternative_solutions.txt results as TSV
"""
import argparse
import multiprocessing
import pandas
import step00


def get_data(file_name: str) -> pandas.DataFrame:
    data = pandas.read_csv(file_name, header=0, sep="\t")
    data.sort_values("SLPP", ascending=False, inplace=True, ignore_index=True)
    data.columns = ["Cellularity", "Ploidy", "SLPP"]
    return data.iloc[0, :]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Sequenza output alternative_solutions.txt file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".txt"), args.input)):
        raise ValueError("INPUT must end with .TXT!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="columns", copy=False).T
        input_data.index = list(map(lambda x: x.split("/")[-2], args.input))
        input_data["Stage"] = pool.map(step00.get_long_sample_type, list(input_data.index))
    print(input_data)

    input_data.to_csv(args.output, sep="\t")
