"""
aggregate_fastqc.py: aggregate FastQC results into a single figure
"""
import argparse
import itertools
import os
import zipfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input ZIP files", type=str, nargs="+")
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--threshold", help="Number of FAIL threshold", type=int, default=0)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".zip"), args.input)):
        raise ValueError("INPUT must end with .zip!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .pdf!!")
    elif args.threshold < 0:
        raise ValueError("Threshold must be non-negative!!")

    args.input.sort(key=step00.sorting)

    for file_name in tqdm.tqdm(args.input):
        with zipfile.ZipFile(file_name) as zip_file:
            zip_file.extract(list(filter(lambda x: x.endswith("summary.txt"), zip_file.namelist()))[0], path=step00.tmpfs)

    raw_data = list()
    for directory in tqdm.tqdm(step00.directory_list(step00.tmpfs)):
        with open(os.path.join(directory, "summary.txt"), "r") as f:
            for line in f.readlines():
                split_line = line.strip().split("\t")
                raw_data.append(split_line)

    data = pandas.DataFrame(raw_data, columns=["Result", "Item", "File"])
    print(data)

    file_set = set(data["File"])
    item_set = set(data["Item"]) - {"Per tile sequence quality", "Per sequence GC content", "Adapter Content"}

    output_data = pandas.DataFrame(numpy.zeros((len(item_set), len(file_set))), index=sorted(item_set), columns=sorted(file_set), dtype=float)
    for f, item in tqdm.tqdm(list(itertools.product(sorted(file_set), sorted(item_set)))):
        result = data.loc[(data["Item"] == item) & (data["File"] == f), "Result"].to_numpy()[0]
        if result == "PASS":
            output_data.loc[item, f] = 1.0
        elif result == "WARN":
            output_data.loc[item, f] = 0.5
        elif result == "FAIL":
            output_data.loc[item, f] = 0.0
        else:
            raise Exception("Something went wrong!!")
    print(output_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(len(file_set) / 5, len(item_set) * 2))

    seaborn.heatmap(data=output_data, vmin=0, center=0.5, vmax=1, cbar=False, square=False, xticklabels=False, cmap="RdYlGn", ax=ax)

    matplotlib.pyplot.xlabel("Files")
    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
