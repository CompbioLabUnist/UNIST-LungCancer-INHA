"""
aggregate_CNVkit_simple.py: Aggregate CNVkit results as simple genome view
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
        return numpy.average(a=a, weights=weights)
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

    size_data = pandas.read_csv(args.size, sep="\t", header=None, names=["chromosome", "length"]).set_index(keys="chromosome", verify_integrity=True)
    print(size_data)

    chromosome_set = set(input_data["chromosome"])
    chromosome_list = list(filter(lambda x: x in chromosome_set, step00.chromosome_full_list))
    print(chromosome_list)

    stage_set = set(map(step00.get_long_sample_type, sample_list))
    stage_list = list(reversed(list(filter(lambda x: x in stage_set, step00.long_sample_type_list))))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, axs = matplotlib.pyplot.subplots(nrows=2 * len(stage_list), ncols=len(chromosome_list), sharex="col", sharey="row", figsize=(len(chromosome_list) * 4, 8 * len(stage_list) * 2), gridspec_kw={"width_ratios": list(map(lambda x: x / step00.big, size_data.loc[chromosome_list, "length"]))})

    for i, chromosome in enumerate(chromosome_list):
        chromosome_data = pandas.DataFrame(data=numpy.ones(shape=(len(sample_list), size_data.loc[chromosome, "length"] // step00.big)), index=sample_list, dtype=float)

        with multiprocessing.Pool(args.cpus) as pool:
            for sample in tqdm.tqdm(sample_list):
                chromosome_data.loc[sample, :] = pool.starmap(get_chromosome_data, [(sample, chromosome, step * step00.big, (step + 1) * step00.big) for step in list(chromosome_data.columns)])

        for j, stage in enumerate(stage_list):
            stage_sample_list = list(filter(lambda x: step00.get_long_sample_type(x) == stage, sample_list))
            proportion = [1.0 for _ in range(chromosome_data.shape[1])]
            if stage_sample_list:
                for k in tqdm.tqdm(range(chromosome_data.shape[1])):
                    proportion[k] = len(list(filter(lambda x: chromosome_data.loc[x, k] >= (1 + args.threshold), stage_sample_list))) / len(stage_sample_list)

            axs[j][i].plot(proportion, color=step00.stage_color_code[stage], linestyle=step00.stage_linestyle[stage], label=stage)

            axs[j][i].set_xticks([])
            axs[j][i].set_ylim(bottom=0, top=1)
            axs[j][i].set_xlabel(chromosome[3:])

            axs[j][i].plot(proportion, color=step00.stage_color_code[stage], linestyle=step00.stage_linestyle[stage], label=stage)

            axs[j][i].set_xticks([])
            axs[j][i].set_ylim(bottom=0, top=1)
            axs[j][i].set_xlabel(chromosome[3:])

            if i == 0:
                axs[j][i].set_ylabel(stage)

        for j, stage in enumerate(stage_list):
            stage_sample_list = list(filter(lambda x: step00.get_long_sample_type(x) == stage, sample_list))
            proportion = [1 for _ in range(chromosome_data.shape[1])]
            if stage_sample_list:
                for k in tqdm.tqdm(range(chromosome_data.shape[1])):
                    proportion[k] = len(list(filter(lambda x: chromosome_data.loc[x, k] <= (1 - args.threshold), stage_sample_list))) / len(stage_sample_list)

            axs[-1 - j][i].plot(proportion, color=step00.stage_color_code[stage], linestyle=step00.stage_linestyle[stage], label=stage)

            axs[-1 - j][i].set_xticks([])
            axs[-1 - j][i].set_ylim(bottom=0, top=1)
            axs[-1 - j][i].set_xlabel(chromosome[3:])
            axs[-1 - j][i].invert_yaxis()

            if i == 0:
                axs[-1 - j][i].set_ylabel(stage)

    matplotlib.pyplot.tight_layout()
    fig.savefig(args.output + ".pdf")
    fig.savefig(args.output + ".png")
    matplotlib.pyplot.close(fig)
