"""
analysis_depth.py: analysis samtools depth
"""
import argparse
import itertools
import typing
import multiprocessing
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import numpy
import pandas
import step00


def get_depth(file_name: str) -> typing.Tuple[str, float, float]:
    """
    get_depth: get depth of coverage from file
    """
    data = pandas.read_csv(file_name, sep="\t", names=["CHROM", "POS", "Depth"], skiprows=1, verbose=True)
    return (step00.get_id(file_name), numpy.mean(data["Depth"]), numpy.std(data["Depth"]) / 2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Samtools depth TSV files", type=str, nargs="+")
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif args.cpus < 0:
        raise ValueError("CPUS must be positive!!")

    args.input.sort(key=step00.sorting)
    IDs = list(map(step00.get_id, args.input))
    print(IDs)

    with multiprocessing.Pool(processes=args.cpus) as pool:
        depth_data = pandas.DataFrame(data=pool.map(get_depth, args.input), columns=["ID", "mean", "std"])
    print(depth_data)

    patient_colors = dict(zip(sorted(set(list(map(step00.get_patient, args.input)))), itertools.cycle(matplotlib.colors.XKCD_COLORS.keys())))
    sample_colors = dict(list(map(lambda x: (x, matplotlib.colors.XKCD_COLORS[patient_colors[step00.get_patient(x)]]), list(map(step00.get_id, args.input)))))
    depth_data["color"] = list(map(lambda x: sample_colors[x], depth_data["ID"]))
    print(depth_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(64, 18))
    matplotlib.pyplot.bar(x=list(depth_data.index), height=depth_data["mean"], color=depth_data["color"], yerr=depth_data["std"])

    matplotlib.pyplot.xticks([])
    matplotlib.pyplot.xlabel("Total {0} samples from {1} Patients".format(len(list(map(step00.get_id, args.input))), len(sorted(set(list(map(step00.get_patient, args.input)))))))
    matplotlib.pyplot.ylabel("Depth")
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
