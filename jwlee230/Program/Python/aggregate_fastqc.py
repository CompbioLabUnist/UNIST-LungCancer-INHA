"""
aggregate_fastqc.py: aggregate FastQC results into a single figure
"""
import argparse
import os
import zipfile
import matplotlib
import matplotlib.pyplot
import pandas
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input ZIP files", type=str, nargs="+")
    parser.add_argument("output", help="Output PNG file", type=str)
    parser.add_argument("--tmpfs", help="Path to temporary directory", type=str, default="/tmpfs")

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".zip"), args.input)):
        raise ValueError("INPUT must end with .zip!!")
    if not args.output.endswith(".png"):
        raise ValueError("Output file must end with .png!!")

    args.input.sort()

    for file_name in args.input:
        with zipfile.ZipFile(file_name) as zip_file:
            zip_file.extract(list(filter(lambda x: x.endswith("summary.txt"), zip_file.namelist()))[0], path=args.tmpfs)

    raw_data = list()
    for directory in step00.directory_list(args.tmpfs):
        with open(os.path.join(directory, "summary.txt"), "r") as f:
            for line in f.readlines():
                split_line = line.strip().split("\t")
                raw_data.append(split_line)

    data = pandas.DataFrame(raw_data, columns=["Result", "Item", "File"])

    item_list = sorted(set(data["Item"]))
    file_list = sorted(set(data["File"]))

    matplotlib.use("Agg")
    matplotlib.rcParams.update({"font.size": 30})

    fig, ax = matplotlib.pyplot.subplots(figsize=(64, 18))

    for i, f in enumerate(file_list):
        for j, item in enumerate(item_list):
            result = list(data.loc[(data["Item"] == item) & (data["File"] == f)]["Result"])[0]

            if result == "PASS":
                matplotlib.pyplot.scatter(i, j, s=20, marker="o", c="g")
            elif result == "WARN":
                matplotlib.pyplot.scatter(i, j, s=20, marker="^", c="y")
            elif result == "FAIL":
                matplotlib.pyplot.scatter(i, j, s=20, marker="X", c="r")
            else:
                raise Exception("Something went wrong!!")

    matplotlib.pyplot.xticks([])
    matplotlib.pyplot.yticks(range(len(item_list)), item_list)

    matplotlib.pyplot.xlabel("Files")
    matplotlib.pyplot.ylabel("Tests")

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
