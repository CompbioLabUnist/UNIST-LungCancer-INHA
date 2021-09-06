"""
aggregate_sequenza_genome.py: Aggregate sequenza results as genome view
"""
import argparse
import multiprocessing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import step00

big = 10 ** 6


def cut_ratio(value: float) -> float:
    if 0 <= value:
        return value
    elif value < 0:
        return 0
    else:
        raise ValueError("Something went wrong!!")


def get_data(file_name: str) -> pandas.DataFrame:
    data = pandas.read_csv(file_name, sep="\t", usecols=["chromosome", "start.pos", "end.pos", "depth.ratio"]).dropna(axis="index")
    data["sample"] = file_name.split("/")[-2]
    data["depth.ratio"] = list(map(cut_ratio, data["depth.ratio"]))
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Sequenza output segments.txt file(s)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("size", help="SIZE file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)
    parser.add_argument("--threshold", help="Threshold for gain/loss", type=float, default=0.2)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".txt"), args.input)):
        raise ValueError("INPUT must end with .TXT!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.threshold < 1):
        raise ValueError("Threshold must be (0, 1)")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    args.input = list(filter(lambda x: step00.get_patient(x.split("/")[-2]) in patients, args.input))

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="index", copy=False, ignore_index=True)
    print(input_data)

    chromosome_list = list(filter(lambda x: x in set(input_data["chromosome"]), step00.chromosome_list))
    sample_list = sorted(set(input_data["sample"]), key=step00.sorting)
    primary_cancer_list = list(filter(lambda x: step00.get_long_sample_type(x) == "Primary", sample_list))
    precancer_list = list(filter(lambda x: step00.get_long_sample_type(x) != "Primary", sample_list))
    print(chromosome_list)
    print(len(sample_list), sample_list)

    size_data = pandas.read_csv(args.size, sep="\t", header=None, names=["chromosome", "length"]).set_index(keys="chromosome", verify_integrity=True)
    print(size_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, axs = matplotlib.pyplot.subplots(nrows=3, ncols=len(chromosome_list), sharex="col", sharey="row", figsize=(len(chromosome_list) * 4, len(sample_list)), gridspec_kw={"height_ratios": [1, len(sample_list) // 5, 1], "width_ratios": list(map(lambda x: x // big, size_data.loc[chromosome_list, "length"]))})
    for i, chromosome in enumerate(chromosome_list):
        print(chromosome)
        chromosome_data = pandas.DataFrame(data=numpy.ones(shape=(len(sample_list), size_data.loc[chromosome, "length"] // big)), index=sample_list, dtype=float)
        for _, row in input_data.loc[(input_data["chromosome"] == chromosome)].iterrows():
            chromosome_data.loc[row["sample"], row["start.pos"] // big:row["end.pos"] // big] = row["depth.ratio"]

        for j in range(chromosome_data.shape[1]):
            axs[0][i].bar(x=j, height=len(list(filter(lambda x: chromosome_data.loc[x, j] >= (1 + args.threshold), primary_cancer_list))) / len(primary_cancer_list), width=1, align="edge", color="tab:red", edgecolor=None, linewidth=0)
            axs[0][i].bar(x=j, height=len(list(filter(lambda x: chromosome_data.loc[x, j] >= (1 + args.threshold), precancer_list))) / len(precancer_list), width=1, align="edge", color="tab:orange", edgecolor=None, linewidth=0)
        axs[0][i].set_ylim(bottom=0, top=1)

        seaborn.heatmap(data=chromosome_data, vmin=0, center=1, vmax=2, cmap="coolwarm", cbar=False, xticklabels=False, yticklabels=True, ax=axs[1][i])
        axs[1][i].set_xlabel(chromosome[3:])

        for j in range(chromosome_data.shape[1]):
            axs[2][i].bar(x=j, height=len(list(filter(lambda x: chromosome_data.loc[x, j] <= (1 - args.threshold), primary_cancer_list))) / len(primary_cancer_list), width=1, align="edge", color="tab:blue", edgecolor=None, linewidth=0)
            axs[2][i].bar(x=j, height=len(list(filter(lambda x: chromosome_data.loc[x, j] <= (1 - args.threshold), precancer_list))) / len(precancer_list), width=1, align="edge", color="tab:cyan", edgecolor=None, linewidth=0)
        axs[2][i].set_ylim(bottom=0, top=1)
        axs[2][i].invert_yaxis()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
