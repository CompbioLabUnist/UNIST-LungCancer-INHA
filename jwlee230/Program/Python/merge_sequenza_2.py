"""
merge_sequenza_2.py: merge Sequenza
"""
import argparse
import multiprocessing
import pandas
import step00


def get_data(file_name: str) -> pandas.DataFrame:
    data = pandas.read_csv(file_name, sep="\t").dropna(axis="index").dropna(axis="index")
    data.columns = ["chromosome", "start", "end"] + list(data.columns)[3:]
    data["Sample"] = file_name.split("/")[-2]
    return data


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
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="index", copy=False, ignore_index=True)
        input_data["Patient"] = pool.map(step00.get_patient, input_data["Sample"])
        input_data["Stage"] = pool.map(step00.get_long_sample_type, input_data["Sample"])
    print(input_data)

    input_data.to_csv(args.output, sep="\t")
