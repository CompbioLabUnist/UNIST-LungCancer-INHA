"""
aggregate_fusioncatcher.py: Aggregate FusionCatcher results as CoMut plot
"""
import argparse
import collections
import multiprocessing
from comut import comut
import matplotlib
import pandas
import step00


def read_maf(filename: str) -> pandas.DataFrame:
    data = pandas.read_csv(filename, sep="\t", low_memory=False)
    data["sample"] = filename.split("/")[-1].split(".")[0]
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="FusionCatcher output .tsv files", type=str, nargs="+")
    parser.add_argument("output", help="Output file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    args.input.sort(key=step00.sorting)
    print(args.input)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    my_comut = comut.CoMut()
    my_comut.samples = sorted(list(map(lambda x: x.split("/")[-1].split(".")[0], args.input)), key=step00.sorting)

    with multiprocessing.Pool(args.cpus) as pool:
        fusioncatcher_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
    print(fusioncatcher_data)

    indicator_data = pandas.DataFrame()
    indicator_data["sample"] = my_comut.samples
    indicator_data["group"] = list(map(lambda x: hash(step00.get_patient(x)), indicator_data["sample"]))
    print(indicator_data)
    my_comut.add_sample_indicators(indicator_data, name="Same patient")
