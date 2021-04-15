"""
aggregate_sequenza_genome.py: Aggregate sequenza results as genome view
"""
import argparse
import multiprocessing
import tarfile
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import step00


def untar(file_name: str) -> pandas.DataFrame:
    with tarfile.open(file_name, "r:gz") as tar:
        txt_file = list(filter(lambda x: x.endswith("_segments.txt"), tar.getnames()))[0]
        print(txt_file)
        tar.extract(txt_file, path=step00.tmpfs)

    data = pandas.read_csv(step00.tmpfs + "/" + txt_file, sep="\t", usecols=["chromosome", "start.pos", "end.pos", "depth.ratio"])
    data.dropna(axis="index", inplace=True)
    data["sample"] = txt_file.split("_")[0]

    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Sequenza tar.gz file(s)", type=str, nargs="+")
    parser.add_argument("size", help="SIZE file", type=str)
    parser.add_argument("centromere", help="Centromere file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tar.gz"), args.input)):
        raise ValueError("INPUT must end with .tar.gz!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    args.input.sort()

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(untar, args.input), axis="index", copy=False, ignore_index=True)
    print(input_data)
    exit()

    colors = list(matplotlib.colors.TABLEAU_COLORS.keys())
    colors_length = len(colors)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, axs = matplotlib.pyplot.subplots(nrows=3, sharex=True, figsize=(18, 64))

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
