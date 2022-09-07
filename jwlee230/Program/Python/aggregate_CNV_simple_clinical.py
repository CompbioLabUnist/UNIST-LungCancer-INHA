"""
aggregate_CNV_genome_clinical.py: Aggregate CNV results as genome view with clinical data
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
watching = ""
threshold = 0.0


def get_chromosome_data(sample: str, chromosome: str, start: int, end: int) -> float:
    tmp_data = input_data.loc[(input_data["Sample"] == sample) & (input_data["chromosome"] == chromosome), :]
    length = end - start + 1

    a = list()
    weights = list()

    get_data = tmp_data.loc[(tmp_data["start"] <= start) & (end <= tmp_data["end"]), :]
    a += list(get_data[watching])
    weights += list(get_data["weight"] * (end - start + 1) / length)

    get_data = tmp_data.loc[(tmp_data["start"] <= start) & (start <= tmp_data["end"]), :]
    a += list(get_data[watching])
    weights += list(get_data["weight"] * (get_data["end"] - start + 1) / length)

    get_data = tmp_data.loc[(tmp_data["start"] <= end) & (end <= tmp_data["end"]), :]
    a += list(get_data[watching])
    weights += list(get_data["weight"] * (end - get_data["start"] + 1) / length)

    get_data = tmp_data.loc[(start <= tmp_data["start"]) & (tmp_data["end"] <= end), :]
    a += list(get_data[watching])
    weights += list(get_data["weight"] * (get_data["end"] - get_data["start"] + 1) / length)

    tmp_start = start
    for index, row in tmp_data.loc[(tmp_data["start"] <= end) & (tmp_data["end"] >= start)].iterrows():
        a.append(1.0)
        weights.append((max(row["start"], start) - tmp_start + 1) / length)
        tmp_start = min(row["end"], end)

    a.append(1.0)
    weights.append((end - tmp_start + 1) / length)

    return 1.0 if ((1 - threshold) < numpy.average(a=a, weights=weights) < (1 + threshold)) else numpy.average(a=a, weights=weights)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="CNV segment.tsv file", type=str)
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("size", help="SIZE file", type=str)
    parser.add_argument("output", help="Output file basename", type=str)
    parser.add_argument("--watching", help="Watching column name", type=str, required=True)
    parser.add_argument("--compare", help="Comparison grouping (type, control, case)", type=str, nargs=3, default=["Recurrence", "NO", "YES"])
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)
    parser.add_argument("--threshold", help="Threshold for gain/loss", type=float, default=0.2)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .TSV!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.threshold < 1):
        raise ValueError("Threshold must be (0, 1)")

    watching = args.watching
    threshold = args.threshold

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

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data = input_data.loc[(input_data["Patient"].isin(patients))]
    print(input_data)

    sample_list = sorted(set(input_data["Sample"]), key=step00.sorting_by_type)
    print(sample_list)

    chromosome_list = list(filter(lambda x: x in set(input_data["chromosome"]), step00.chromosome_list))
    print(chromosome_list)

    control_sample_list = list(filter(lambda x: step00.get_patient(x) in control_patients, sample_list))
    control_primary_list = list(filter(lambda x: step00.get_long_sample_type(x) == "Primary", control_sample_list))
    control_precancer_list = list(filter(lambda x: step00.get_long_sample_type(x) != "Primary", control_sample_list))

    case_sample_list = list(filter(lambda x: step00.get_patient(x) in case_patients, sample_list))
    case_primary_list = list(filter(lambda x: step00.get_long_sample_type(x) == "Primary", case_sample_list))
    case_precancer_list = list(filter(lambda x: step00.get_long_sample_type(x) != "Primary", case_sample_list))

    size_data = pandas.read_csv(args.size, sep="\t", header=None, names=["chromosome", "length"]).set_index(keys="chromosome", verify_integrity=True)
    print(size_data)

    stage_set = set(map(step00.get_long_sample_type, sample_list))
    stage_list = list(filter(lambda x: x in stage_set, reversed(step00.long_sample_type_list)))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, axs = matplotlib.pyplot.subplots(nrows=2 * len(stage_list), ncols=len(chromosome_list), sharex="col", sharey="row", figsize=(len(chromosome_list) * 4, 8 * len(stage_list) * 2), gridspec_kw={"width_ratios": list(map(lambda x: x / step00.big, size_data.loc[chromosome_list, "length"]))})

    for i, chromosome in enumerate(chromosome_list):
        chromosome_data = pandas.DataFrame(data=numpy.ones(shape=(len(control_sample_list + case_sample_list), size_data.loc[chromosome, "length"] // step00.big)), index=control_sample_list + case_sample_list, dtype=float)

        with multiprocessing.Pool(args.cpus) as pool:
            for sample in tqdm.tqdm(sample_list):
                chromosome_data.loc[sample, :] = pool.starmap(get_chromosome_data, [(sample, chromosome, step * step00.big, (step + 1) * step00.big) for step in list(chromosome_data.columns)])

        for j, stage in enumerate(stage_list):
            control_stage_sample_list = list(filter(lambda x: step00.get_long_sample_type(x) == stage, control_sample_list))
            case_stage_sample_list = list(filter(lambda x: step00.get_long_sample_type(x) == stage, case_sample_list))

            control_proportion = [1.0 for _ in range(chromosome_data.shape[1])]
            case_proportion = [1.0 for _ in range(chromosome_data.shape[1])]
            for k in tqdm.tqdm(range(chromosome_data.shape[1])):
                if control_stage_sample_list:
                    control_proportion[k] = len(list(filter(lambda x: chromosome_data.loc[x, k] >= (1 + args.threshold), control_stage_sample_list))) / len(control_stage_sample_list)
                if case_stage_sample_list:
                    case_proportion[k] = len(list(filter(lambda x: chromosome_data.loc[x, k] >= (1 + args.threshold), case_stage_sample_list))) / len(case_stage_sample_list)

            axs[j][i].plot(control_proportion, color=step00.stage_color_code[stage], linestyle="--", label=args.compare[1])
            axs[j][i].plot(case_proportion, color=step00.stage_color_code[stage], linestyle="-", label=args.compare[2])

            axs[j][i].set_xticks([])
            axs[j][i].set_ylim(bottom=0, top=1)
            axs[j][i].set_xlabel(chromosome[3:])

            if i == 0:
                axs[j][i].set_ylabel(stage)
                axs[j][i].legend(title=args.compare[0], loc="upper center")

        for j, stage in enumerate(stage_list):
            control_stage_sample_list = list(filter(lambda x: step00.get_long_sample_type(x) == stage, control_sample_list))
            case_stage_sample_list = list(filter(lambda x: step00.get_long_sample_type(x) == stage, case_sample_list))

            control_proportion = [1.0 for _ in range(chromosome_data.shape[1])]
            case_proportion = [1.0 for _ in range(chromosome_data.shape[1])]
            for k in tqdm.tqdm(range(chromosome_data.shape[1])):
                if control_stage_sample_list:
                    control_proportion[k] = len(list(filter(lambda x: chromosome_data.loc[x, k] <= (1 - args.threshold), control_stage_sample_list))) / len(control_stage_sample_list)
                if case_stage_sample_list:
                    case_proportion[k] = len(list(filter(lambda x: chromosome_data.loc[x, k] <= (1 - args.threshold), case_stage_sample_list))) / len(case_stage_sample_list)

            axs[-1 - j][i].plot(control_proportion, color=step00.stage_color_code[stage], linestyle="--", label=args.compare[1])
            axs[-1 - j][i].plot(case_proportion, color=step00.stage_color_code[stage], linestyle="-", label=args.compare[2])

            axs[-1 - j][i].set_xticks([])
            axs[-1 - j][i].set_ylim(bottom=0, top=1)
            axs[-1 - j][i].set_xlabel(chromosome[3:])
            axs[-1 - j][i].invert_yaxis()

            if i == 0:
                axs[-1 - j][i].set_ylabel(stage)
                axs[-1 - j][i].legend(title=args.compare[0], loc="lower center")

    matplotlib.pyplot.tight_layout()
    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
