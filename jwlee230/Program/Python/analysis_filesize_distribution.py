"""
analysis_filesize_distribution.py: analysis file-size distribution for QC
"""
import argparse
import itertools
import os
import multiprocessing
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import numpy
import pandas
import step00


def get_size(file):
    return (file, os.path.getsize(file) / (1024 ** 3))


def get_id(file):
    return step00.get_id(file).split("_")[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input FASTQ.gz file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .pdf!!")
    elif args.cpus < 0:
        raise ValueError("CPUS must be positive!!")

    args.input.sort()

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.DataFrame(data=pool.map(get_size, args.input), columns=["File", "Size"])
        input_data["ID"] = pool.map(get_id, input_data["File"])
    print(input_data)

    patient_colors = dict(zip(sorted(set(list(map(step00.get_patient, input_data["ID"])))), itertools.cycle(matplotlib.colors.XKCD_COLORS.keys())))
    sample_colors = dict(list(map(lambda x: (x, matplotlib.colors.XKCD_COLORS[patient_colors[step00.get_patient(x)]]), input_data["ID"])))
    input_data["color"] = list(map(lambda x: sample_colors[x], input_data["ID"]))
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(64, 18))

    matplotlib.pyplot.bar(x=list(input_data.index), height=input_data["Size"], color=input_data["color"])
    matplotlib.pyplot.axhline(numpy.mean(input_data["Size"]), color="k", linestyle="--")
    matplotlib.pyplot.text(0, numpy.mean(input_data["Size"]), f"Mean: {numpy.mean(input_data['Size']):.1f} GB", color="k", horizontalalignment="left", verticalalignment="baseline")
    matplotlib.pyplot.axhline(numpy.median(input_data["Size"]), color="k", linestyle=":")
    matplotlib.pyplot.text(len(input_data), numpy.median(input_data["Size"]), f"Median: {numpy.median(input_data['Size']):.1f} GB", color="k", horizontalalignment="right", verticalalignment="baseline")

    matplotlib.pyplot.xticks([])
    matplotlib.pyplot.xlabel("Total {0} paired samples from {1} Patients".format(len(list(map(step00.get_id, args.input))) / 2, len(sorted(set(list(map(step00.get_patient, args.input)))))))
    matplotlib.pyplot.ylabel("Size (GB)")
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
