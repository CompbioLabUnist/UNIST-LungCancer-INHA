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

    CNV_dict = dict()
    CNV_length = dict()
    for c in step00.chromosome_list:
        CNV_dict[c] = 0.0
        CNV_length[c] = 0

    for i, row in data.iterrows():
        length = row["end.pos"] - row["start.pos"] + 1
        CNV_dict[row["chromosome"]] += row["depth.ratio"] * length
        CNV_length[row["chromosome"]] += length

    for c in step00.chromosome_list:
        CNV_dict[c] /= CNV_length[c]

    return pandas.DataFrame.from_dict(CNV_dict, orient="index", columns=[txt_file.split("_")[0]])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Sequenza tar.gz file(s)", type=str, nargs="+")
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
        input_data = pandas.concat(objs=pool.map(untar, args.input), axis="columns", copy=False, verify_integrity=True)
    input_data = input_data[sorted(input_data.columns, key=step00.sorting)]
    print(input_data)

    colors = list(matplotlib.colors.TABLEAU_COLORS.keys())
    colors_length = len(colors)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(64, 18))

    seaborn.heatmap(data=input_data, cmap="vlag", vmin=0, center=1, vmax=2, cbar=False, xticklabels=False, yticklabels=True, linewidth=0.1, linecolor="white", ax=ax)
    matplotlib.pyplot.text(0.5, input_data.shape[0] - 0.5, "Same patient", fontsize="xx-small", horizontalalignment="left", verticalalignment="bottom", color="black")

    patient_IDs = list(map(step00.get_patient, list(input_data.columns)))
    for i, ID in enumerate(sorted(set(patient_IDs))):
        xs = list(map(lambda x: x + 0.5, step00.list_first_last(patient_IDs, ID)))
        matplotlib.pyplot.plot(numpy.arange(xs[0], xs[1] + 1), [input_data.shape[0] - 0.2 for _ in numpy.arange(xs[0], xs[1] + 1)], "o-", color=matplotlib.colors.TABLEAU_COLORS[colors[i % colors_length]], markersize=7)

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
