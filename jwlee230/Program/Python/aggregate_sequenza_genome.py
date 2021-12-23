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
import tqdm
import step00

watching = "depth.ratio"


def get_data(file_name: str) -> pandas.DataFrame:
    data = pandas.read_csv(file_name, sep="\t", usecols=["chromosome", "start.pos", "end.pos", watching]).dropna(axis="index")
    data["sample"] = file_name.split("/")[-2]
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Sequenza output segments.txt file(s)", type=str, nargs="+")
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

    if list(filter(lambda x: not x.endswith(".txt"), args.input)):
        raise ValueError("INPUT must end with .TXT!!")
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

    args.input = list(filter(lambda x: step00.get_patient(x.split("/")[-2]) in patients, args.input))

    if args.patient:
        args.input.sort(key=lambda x: step00.sorting(x.split("/")[-2]))
    elif args.type:
        args.input.sort(key=lambda x: step00.sorting_by_type(x.split("/")[-2]))
    else:
        raise Exception("Something went wrong!!")

    sample_list = list(map(lambda x: x.split("/")[-2], args.input))

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="index", copy=False, ignore_index=True, verify_integrity=True)
    print(input_data)

    chromosome_list = list(filter(lambda x: x in set(input_data["chromosome"]), step00.chromosome_full_list))
    primary_cancer_list = list(filter(lambda x: step00.get_long_sample_type(x) == "Primary", sample_list))
    precancer_list = list(filter(lambda x: step00.get_long_sample_type(x) != "Primary", sample_list))
    print(chromosome_list)

    size_data = pandas.read_csv(args.size, sep="\t", header=None, names=["chromosome", "length"]).set_index(keys="chromosome", verify_integrity=True)
    print(size_data)

    stage_set = set(map(step00.get_long_sample_type, sample_list))
    stage_list = list(filter(lambda x: x in stage_set, step00.long_sample_type_list))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, axs = matplotlib.pyplot.subplots(nrows=3, ncols=len(chromosome_list), sharex="col", sharey="row", figsize=(len(chromosome_list) * 4, len(sample_list) + 16), gridspec_kw={"height_ratios": [8, len(sample_list), 8], "width_ratios": list(map(lambda x: x / step00.big, size_data.loc[chromosome_list, "length"]))})

    for i, chromosome in enumerate(chromosome_list):
        chromosome_data = pandas.DataFrame(data=numpy.ones(shape=(len(sample_list), size_data.loc[chromosome, "length"] // step00.big)), index=sample_list, dtype=float)

        for index, row in tqdm.tqdm(input_data.loc[(input_data["chromosome"] == chromosome)].iterrows()):
            chromosome_data.loc[row["sample"], row["start.pos"] // step00.big:row["end.pos"] // step00.big] = row[watching]

        for stage in stage_list:
            stage_sample_list = list(filter(lambda x: step00.get_long_sample_type(x) == stage, sample_list))
            proportion = [1 for _ in range(chromosome_data.shape[1])]

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
            proportion = [1 for _ in range(chromosome_data.shape[1])]

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
