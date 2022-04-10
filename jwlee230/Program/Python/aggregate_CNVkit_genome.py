"""
aggregate_CNVkit_genome.py: Aggregate CNVkit results as genome view
"""
import argparse
import multiprocessing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import tqdm
import step00

input_data = pandas.DataFrame()
watching = "log2"


def get_data(file_name: str) -> pandas.DataFrame:
    data = pandas.read_csv(file_name, sep="\t")
    data["ID"] = step00.get_id(file_name)
    return data


def get_chromosome_data(sample: str, chromosome: str, start: int, end: int) -> float:
    tmp_data = input_data.loc[(input_data["ID"] == sample) & (input_data["chromosome"] == chromosome), :]

    a = list()
    weights = list()

    get_data = tmp_data.loc[(tmp_data["start"] <= start) & (end <= tmp_data["end"]), :]
    a += list(get_data["exp"])
    weights += list(get_data["weight"] * (end - start + 1))

    get_data = tmp_data.loc[(tmp_data["start"] <= start) & (start <= tmp_data["end"]), :]
    a += list(get_data["exp"])
    weights += list(get_data["weight"] * (get_data["end"] - start + 1))

    get_data = tmp_data.loc[(tmp_data["start"] <= end) & (end <= tmp_data["end"]), :]
    a += list(get_data["exp"])
    weights += list(get_data["weight"] * (end - get_data["start"] + 1))

    get_data = tmp_data.loc[(start <= tmp_data["start"]) & (tmp_data["end"] <= end), :]
    a += list(get_data["exp"])
    weights += list(get_data["weight"] * (get_data["end"] - get_data["start"] + 1))

    if a and weights:
        return 1.0 if (0.8 < numpy.average(a=a, weights=weights) < 1.2) else numpy.average(a=a, weights=weights)
    else:
        return 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PureCN output CNR file(s)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("size", help="SIZE file", type=str)
    parser.add_argument("output", help="Output file basename", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)
    parser.add_argument("--threshold", help="Threshold for gain/loss", type=float, default=0.2)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_sorting = parser.add_mutually_exclusive_group(required=True)
    group_sorting.add_argument("--patient", help="Sorting by patient first", action="store_true", default=False)
    group_sorting.add_argument("--type", help="Sorting by type first", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".cnr"), args.input)):
        raise ValueError("INPUT must end with .CNR!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
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

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))

    if args.patient:
        args.input.sort(key=step00.sorting)
    elif args.type:
        args.input.sort(key=step00.sorting_by_type)
    else:
        raise Exception("Something went wrong!!")

    sample_list = list(map(step00.get_id, args.input))
    print(len(sample_list), sample_list)

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="index", copy=False, ignore_index=True, verify_integrity=True)
    input_data["exp"] = numpy.power(2, input_data[watching])
    print(input_data)

    chromosome_set = set(input_data["chromosome"])
    chromosome_list = list(filter(lambda x: x in chromosome_set, step00.chromosome_list))
    print(chromosome_list)

    size_data = pandas.read_csv(args.size, sep="\t", header=None, names=["chromosome", "length"]).set_index(keys="chromosome", verify_integrity=True)
    print(size_data)

    stage_set = set(map(step00.get_long_sample_type, args.input))
    stage_list = list(filter(lambda x: x in stage_set, step00.long_sample_type_list))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, axs = matplotlib.pyplot.subplots(nrows=3, ncols=len(chromosome_list), sharex="col", sharey="row", figsize=(len(chromosome_list) * 4, len(sample_list) + 16), gridspec_kw={"height_ratios": [8, len(sample_list), 8], "width_ratios": list(map(lambda x: x / step00.big, size_data.loc[chromosome_list, "length"]))})

    for i, chromosome in enumerate(chromosome_list):
        chromosome_data = pandas.DataFrame(data=numpy.ones(shape=(len(sample_list), size_data.loc[chromosome, "length"] // step00.big)), index=sample_list, dtype=float)

        with multiprocessing.Pool(args.cpus) as pool:
            for sample in tqdm.tqdm(sample_list):
                chromosome_data.loc[sample, :] = pool.starmap(get_chromosome_data, [(sample, chromosome, step * step00.big, (step + 1) * step00.big) for step in list(chromosome_data.columns)])

        for stage in stage_list:
            stage_sample_list = list(filter(lambda x: step00.get_long_sample_type(x) == stage, sample_list))
            proportion = [1.0 for _ in range(chromosome_data.shape[1])]
            for j in tqdm.tqdm(range(chromosome_data.shape[1])):
                proportion[j] = len(list(filter(lambda x: chromosome_data.loc[x, j] >= (1 + args.threshold), stage_sample_list))) / len(stage_sample_list)
            axs[0][i].plot(proportion, color=step00.stage_color_code[stage], linestyle=step00.stage_linestyle[stage], label=stage)

        axs[0][i].set_ylim(bottom=0, top=1)
        axs[0][i].set_xticks([])
        axs[0][i].set_xlabel(chromosome[3:])

        if i == 0:
            axs[0][i].set_ylabel("Proportion")
            axs[0][i].legend(loc="upper center")

        seaborn.heatmap(data=chromosome_data, vmin=0, center=1, vmax=2, cmap="coolwarm", cbar=False, xticklabels=False, yticklabels=True, ax=axs[1][i])
        axs[1][i].set_xlabel(chromosome[3:])

        for stage in stage_list:
            stage_sample_list = list(filter(lambda x: step00.get_long_sample_type(x) == stage, sample_list))
            proportion = [1.0 for _ in range(chromosome_data.shape[1])]
            for j in tqdm.tqdm(range(chromosome_data.shape[1])):
                proportion[j] = len(list(filter(lambda x: chromosome_data.loc[x, j] <= (1 - args.threshold), stage_sample_list))) / len(stage_sample_list)
            axs[2][i].plot(proportion, color=step00.stage_color_code[stage], linestyle=step00.stage_linestyle[stage], label=stage)

        axs[2][i].set_ylim(bottom=0, top=1)
        axs[2][i].invert_yaxis()
        axs[2][i].set_xticks([])
        axs[2][i].set_xlabel(chromosome[3:])

        if i == 0:
            axs[2][i].set_ylabel("Proportion")
            axs[2][i].legend(loc="lower center")

    matplotlib.pyplot.tight_layout()
    fig.savefig(args.output + ".pdf")
    fig.savefig(args.output + ".png")
    matplotlib.pyplot.close(fig)
