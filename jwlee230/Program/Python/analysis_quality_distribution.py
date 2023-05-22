"""
analysis_quality_distribution.py: Analysis of quality distribution by Picard
"""
import argparse
import itertools
import math
import typing
import multiprocessing
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import numpy
import pandas
import step00


def get_quality(file_name: str) -> typing.Tuple[str, float, float, int]:
    """
    get_quality: get depth of coverage from file
    """
    data = pandas.read_csv(file_name, sep="\t", comment="#", dtype=int)
    mean = numpy.average(data["QUALITY"], weights=data["COUNT_OF_Q"])
    std = math.sqrt(numpy.average((data["QUALITY"] - mean) ** 2, weights=data["COUNT_OF_Q"]))
    return (step00.get_id(file_name), mean, std, sum(data["COUNT_OF_Q"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Samtools depth TSV files (not necessarily TSV)", type=str, nargs="+")
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .pdf!!")
    elif args.cpus < 0:
        raise ValueError("CPUS must be positive!!")

    args.input.sort(key=step00.sorting)
    IDs = list(map(step00.get_id, args.input))
    print(IDs)

    with multiprocessing.Pool(processes=args.cpus) as pool:
        quality_data = pandas.DataFrame(data=pool.map(get_quality, args.input), columns=["ID", "mean", "std", "num"])
    print(quality_data)

    mean_value = sum(list(map(lambda x: x[1] * x[3], quality_data.itertuples(index=False, name=None)))) / sum(quality_data["num"])

    patient_colors = dict(zip(sorted(set(list(map(step00.get_patient, args.input)))), itertools.cycle(matplotlib.colors.XKCD_COLORS.keys())))
    sample_colors = dict(list(map(lambda x: (x, matplotlib.colors.XKCD_COLORS[patient_colors[step00.get_patient(x)]]), list(map(step00.get_id, args.input)))))
    quality_data["color"] = list(map(lambda x: sample_colors[x], quality_data["ID"]))
    print(quality_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(64, 18))

    matplotlib.pyplot.bar(x=list(quality_data.index), height=quality_data["mean"], color=quality_data["color"], yerr=quality_data["std"])

    matplotlib.pyplot.axhline(mean_value, color="k", linestyle="--")
    matplotlib.pyplot.text(0, mean_value, f"Mean: {mean_value:.1f}", color="k", horizontalalignment="left", verticalalignment="baseline")

    matplotlib.pyplot.xticks([])
    matplotlib.pyplot.xlabel("Total {0} samples from {1} Patients".format(len(list(map(step00.get_id, args.input))), len(sorted(set(list(map(step00.get_patient, args.input)))))))
    matplotlib.pyplot.ylabel("Quality")
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
