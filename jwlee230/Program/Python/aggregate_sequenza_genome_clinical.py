"""
aggregate_sequenza_genome_clinical.py: Aggregate sequenza results as genome view with clinical data
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
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--compare", help="Comparison grouping (type, control, case)", type=str, nargs=3, default=["Recurrence", "NO", "YES"])
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
        control_patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC") & (clinical_data[args.compare[0]] == args.compare[1])].index)
        case_patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC") & (clinical_data[args.compare[0]] == args.compare[2])].index)
    elif args.ADC:
        control_patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC") & (clinical_data[args.compare[0]] == args.compare[1])].index)
        case_patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC") & (clinical_data[args.compare[0]] == args.compare[2])].index)
    else:
        raise Exception("Something went wrong!!")
    patients = set(list(control_patients) + list(case_patients))
    print(sorted(control_patients))
    print(sorted(case_patients))

    args.input = list(filter(lambda x: step00.get_patient(x.split("/")[-2]) in control_patients, args.input)) + list(filter(lambda x: step00.get_patient(x.split("/")[-2]) in case_patients, args.input))

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="index", copy=False, ignore_index=True, verify_integrity=True)
    print(input_data)

    chromosome_list = list(filter(lambda x: x in set(input_data["chromosome"]), step00.chromosome_list))
    print(chromosome_list)

    control_sample_list = sorted(list(filter(lambda x: step00.get_patient(x) in control_patients, set(input_data["sample"]))), key=step00.sorting_by_type)
    control_primary_list = list(filter(lambda x: step00.get_long_sample_type(x) == "Primary", control_sample_list))
    control_precancer_list = list(filter(lambda x: step00.get_long_sample_type(x) != "Primary", control_sample_list))

    case_sample_list = sorted(list(filter(lambda x: step00.get_patient(x) in case_patients, set(input_data["sample"]))), key=step00.sorting_by_type)
    case_primary_list = list(filter(lambda x: step00.get_long_sample_type(x) == "Primary", case_sample_list))
    case_precancer_list = list(filter(lambda x: step00.get_long_sample_type(x) != "Primary", case_sample_list))

    size_data = pandas.read_csv(args.size, sep="\t", header=None, names=["chromosome", "length"]).set_index(keys="chromosome", verify_integrity=True)
    print(size_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, axs = matplotlib.pyplot.subplots(nrows=6, ncols=len(chromosome_list), sharex="col", sharey="row", figsize=(len(chromosome_list) * 4, len(control_sample_list + case_sample_list) + 32), gridspec_kw={"height_ratios": [8, 8, len(control_sample_list), len(case_sample_list), 8, 8], "width_ratios": list(map(lambda x: x / step00.big, size_data.loc[chromosome_list, "length"]))})

    for i, chromosome in enumerate(chromosome_list):
        chromosome_data = pandas.DataFrame(data=numpy.ones(shape=(len(control_sample_list + case_sample_list), size_data.loc[chromosome, "length"] // step00.big)), index=control_sample_list + case_sample_list, dtype=float)

        for _, row in tqdm.tqdm(input_data.loc[(input_data["chromosome"] == chromosome)].iterrows()):
            chromosome_data.loc[row["sample"], row["start.pos"] // step00.big:row["end.pos"] // step00.big] = row[watching]

        control_primary_proportion = list()
        control_precancer_proportion = list()
        case_primary_proportion = list()
        case_precancer_proportion = list()
        for j in tqdm.tqdm(range(chromosome_data.shape[1])):
            control_primary_proportion.append(len(list(filter(lambda x: chromosome_data.loc[x, j] >= (1 + args.threshold), control_primary_list))) / len(control_primary_list))
            control_precancer_proportion.append(len(list(filter(lambda x: chromosome_data.loc[x, j] >= (1 + args.threshold), control_precancer_list))) / len(control_precancer_list))
            case_primary_proportion.append(len(list(filter(lambda x: chromosome_data.loc[x, j] >= (1 + args.threshold), case_primary_list))) / len(case_primary_list))
            case_precancer_proportion.append(len(list(filter(lambda x: chromosome_data.loc[x, j] >= (1 + args.threshold), case_precancer_list))) / len(case_precancer_list))

        axs[0][i].plot(range(chromosome_data.shape[1]), control_primary_proportion, color="red", linestyle="--", label=args.compare[1])
        axs[0][i].plot(range(chromosome_data.shape[1]), case_primary_proportion, color="red", linestyle="-", label=args.compare[2])
        axs[0][i].set_ylim(bottom=0, top=1)
        axs[0][i].set_xlabel(chromosome[3:])
        if i == 0:
            axs[0][i].set_ylabel("Primary")
            axs[0][i].legend(title=args.compare[0], loc="upper center")

        axs[1][i].plot(range(chromosome_data.shape[1]), control_precancer_proportion, color="lightsalmon", linestyle="--", label=args.compare[1])
        axs[1][i].plot(range(chromosome_data.shape[1]), case_precancer_proportion, color="lightsalmon", linestyle="-", label=args.compare[2])
        axs[1][i].set_ylim(bottom=0, top=1)
        axs[1][i].set_xlabel(chromosome[3:])
        if i == 0:
            axs[1][i].set_ylabel("Precancer")
            axs[1][i].legend(title=args.compare[0], loc="upper center")

        seaborn.heatmap(data=chromosome_data.loc[control_sample_list, :], vmin=0, center=1, vmax=2, cmap="coolwarm", cbar=False, xticklabels=False, yticklabels=True, ax=axs[2][i])
        axs[2][i].set_xlabel(chromosome[3:])
        if i == 0:
            axs[2][i].set_ylabel("{0} - {1}".format(args.compare[0], args.compare[1]))

        seaborn.heatmap(data=chromosome_data.loc[case_sample_list, :], vmin=0, center=1, vmax=2, cmap="coolwarm", cbar=False, xticklabels=False, yticklabels=True, ax=axs[3][i])
        axs[3][i].set_xlabel(chromosome[3:])
        if i == 0:
            axs[3][i].set_ylabel("{0} - {1}".format(args.compare[0], args.compare[2]))

        control_primary_proportion = list()
        control_precancer_proportion = list()
        case_primary_proportion = list()
        case_precancer_proportion = list()
        for j in tqdm.tqdm(range(chromosome_data.shape[1])):
            control_primary_proportion.append(len(list(filter(lambda x: chromosome_data.loc[x, j] <= (1 - args.threshold), control_primary_list))) / len(control_primary_list))
            control_precancer_proportion.append(len(list(filter(lambda x: chromosome_data.loc[x, j] <= (1 - args.threshold), control_precancer_list))) / len(control_precancer_list))
            case_primary_proportion.append(len(list(filter(lambda x: chromosome_data.loc[x, j] <= (1 - args.threshold), case_primary_list))) / len(case_primary_list))
            case_precancer_proportion.append(len(list(filter(lambda x: chromosome_data.loc[x, j] <= (1 - args.threshold), case_precancer_list))) / len(case_precancer_list))

        axs[4][i].plot(range(chromosome_data.shape[1]), control_precancer_proportion, color="cyan", linestyle="--", label=args.compare[1])
        axs[4][i].plot(range(chromosome_data.shape[1]), case_precancer_proportion, color="cyan", linestyle="-", label=args.compare[2])
        axs[4][i].set_ylim(bottom=0, top=1)
        axs[4][i].invert_yaxis()
        axs[4][i].set_xticks([])
        axs[4][i].set_xlabel(chromosome[3:])
        if i == 0:
            axs[4][i].set_ylabel("Precancer")
            axs[4][i].legend(title=args.compare[0], loc="lower center")

        axs[5][i].plot(range(chromosome_data.shape[1]), control_primary_proportion, color="navy", linestyle="--", label=args.compare[1])
        axs[5][i].plot(range(chromosome_data.shape[1]), case_primary_proportion, color="navy", linestyle="-", label=args.compare[2])
        axs[5][i].set_ylim(bottom=0, top=1)
        axs[5][i].invert_yaxis()
        axs[5][i].set_xticks([])
        axs[5][i].set_xlabel(chromosome[3:])
        if i == 0:
            axs[5][i].set_ylabel("Primary")
            axs[5][i].legend(title=args.compare[0], loc="lower center")

    matplotlib.pyplot.tight_layout()
    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
