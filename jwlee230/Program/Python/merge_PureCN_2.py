"""
merge_PureCN.py: merge PureCN csv results as TSV
"""
import argparse
import multiprocessing
import pandas
import step00


def get_data(file_name: str) -> pandas.DataFrame:
    data = pandas.read_csv(file_name, sep="\t")
    data["chrom"] = list(map(lambda x: "chrX" if (str(x) == "23") else ("chr" + str(x)), data["chrom"]))
    data.columns = ["Sample", "chromosome", "start", "end"] + list(data.columns)[4:]
    data["Sample"] = step00.get_id(file_name)
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PureCN seg.TSV file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="index", copy=False, ignore_index=True)
        print(input_data)
        input_data["Patient"] = pool.map(step00.get_patient, input_data["Sample"])
        input_data["Stage"] = pool.map(step00.get_long_sample_type, input_data["Sample"])
    input_data["weight"] = 1.0
    print(input_data)

    input_data.to_csv(args.output, sep="\t")
