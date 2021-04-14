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
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--threshold", help="Number of FAIL threshold", type=int, default=0)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".zip"), args.input)):
        raise ValueError("INPUT must end with .zip!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .pdf!!")
    elif args.threshold < 0:
        raise ValueError("Threshold must be positive!!")

    args.input.sort()

    for file_name in args.input:
        with zipfile.ZipFile(file_name) as zip_file:
            zip_file.extract(list(filter(lambda x: x.endswith("summary.txt"), zip_file.namelist()))[0], path=step00.tmpfs)

    raw_data = list()
    for directory in step00.directory_list(step00.tmpfs):
        with open(os.path.join(directory, "summary.txt"), "r") as f:
            for line in f.readlines():
                split_line = line.strip().split("\t")
                raw_data.append(split_line)

    data = pandas.DataFrame(raw_data, columns=["Result", "Item", "File"])

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(64, 18))

    for i, f in enumerate(sorted(set(data["File"]))):
        FAIL_count = 0
        for j, item in enumerate(sorted(set(data["Item"]))):
            result = data.loc[(data["Item"] == item) & (data["File"] == f), "Result"].to_numpy()[0]

            if result == "PASS":
                matplotlib.pyplot.scatter(i, j, s=30, marker="o", c="tab:green")
            elif result == "WARN":
                matplotlib.pyplot.scatter(i, j, s=30, marker="^", c="tab:olive")
            elif result == "FAIL":
                matplotlib.pyplot.scatter(i, j, s=30, marker="X", c="tab:red")
                FAIL_count += 1
            else:
                raise Exception("Something went wrong!!")
        if FAIL_count > args.threshold:
            print(f, FAIL_count)

    matplotlib.pyplot.xticks([])
    matplotlib.pyplot.yticks(range(len(set(data["Item"]))), sorted(set(data["Item"])), fontsize="xx-small")
    matplotlib.pyplot.xlabel("Files")
    matplotlib.pyplot.ylabel("Tests")

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
