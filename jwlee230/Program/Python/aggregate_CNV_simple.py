"""
aggregate_CNV_simple.py: Aggregate CNV results as simple view
"""
import argparse
import typing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="CNV segment.tsv file", type=str)
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("size", help="SIZE file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--watching", help="Watching column name", type=str, required=True)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)
    parser.add_argument("--threshold", help="Threshold for gain/loss", type=float, default=0.2)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
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

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data = input_data.loc[(input_data["Patient"].isin(patients))]
    print(input_data)

    sample_list = sorted(set(input_data["Sample"]), key=step00.sorting)
    stage_dict: typing.Dict[str, typing.List[str]] = {x: [] for x in step00.long_sample_type_list}
    for sample in tqdm.tqdm(sample_list):
        stage_dict[step00.get_long_sample_type(sample)].append(sample)
    input_list = list(filter(None, list(stage_dict.values())))

    print(list(map(len, input_list)))
    print(input_list)

    chromosome_list = list(filter(lambda x: x in set(input_data["chromosome"]), step00.chromosome_list))
    primary_cancer_list = list(filter(lambda x: step00.get_long_sample_type(x) == "Primary", sample_list))
    precancer_list = list(filter(lambda x: step00.get_long_sample_type(x) != "Primary", sample_list))
    print(chromosome_list)

    size_data = pandas.read_csv(args.size, sep="\t", header=None, names=["chromosome", "length"]).set_index(keys="chromosome", verify_integrity=True)
    print(size_data)

    stage_set = set(map(step00.get_long_sample_type, sample_list))
    stage_list = list(filter(lambda x: x in stage_set, reversed(step00.long_sample_type_list)))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, axs = matplotlib.pyplot.subplots(nrows=2 * len(stage_list), ncols=len(chromosome_list), sharex="col", sharey="row", figsize=(len(chromosome_list) * 4, 8 * len(stage_list) * 2), gridspec_kw={"width_ratios": list(map(lambda x: x / step00.big, size_data.loc[chromosome_list, "length"]))})

    for i, chromosome in enumerate(chromosome_list):
        chromosome_data = pandas.DataFrame(data=numpy.ones(shape=(len(sample_list), size_data.loc[chromosome, "length"] // step00.big)), index=sample_list, dtype=float)

        for index, row in tqdm.tqdm(input_data.loc[(input_data["chromosome"] == chromosome)].iterrows()):
            chromosome_data.loc[row["Sample"], row["start"] // step00.big:row["end"] // step00.big] = 1.0 if (0.8 < row[args.watching] < 1.2) else row[args.watching]

        for j, stage in enumerate(stage_list):
            stage_sample_list = list(filter(lambda x: step00.get_long_sample_type(x) == stage, sample_list))
            proportion = [1.0 for _ in range(chromosome_data.shape[1])]
            if stage_sample_list:
                for k in tqdm.tqdm(range(chromosome_data.shape[1])):
                    proportion[k] = len(list(filter(lambda x: chromosome_data.loc[x, k] >= (1 + args.threshold), stage_sample_list))) / len(stage_sample_list)

            axs[j][i].plot(proportion, color=step00.stage_color_code[stage], linestyle=step00.stage_linestyle[stage], label=stage)
            axs[j][i].fill_between(range(len(proportion)), proportion, color=step00.stage_color_code[stage])

            axs[j][i].set_xticks([])
            axs[j][i].set_ylim(bottom=0, top=1)
            axs[j][i].set_xlabel(chromosome[3:])

            if i == 0:
                axs[j][i].set_ylabel(stage)

        for j, stage in enumerate(stage_list):
            stage_sample_list = list(filter(lambda x: step00.get_long_sample_type(x) == stage, sample_list))
            proportion = [1.0 for _ in range(chromosome_data.shape[1])]
            if stage_sample_list:
                for k in tqdm.tqdm(range(chromosome_data.shape[1])):
                    proportion[k] = len(list(filter(lambda x: chromosome_data.loc[x, k] <= (1 - args.threshold), stage_sample_list))) / len(stage_sample_list)

            axs[-1 - j][i].plot(proportion, color=step00.stage_color_code[stage], linestyle=step00.stage_linestyle[stage], label=stage)
            axs[-1 - j][i].fill_between(range(len(proportion)), proportion, color=step00.stage_color_code[stage])

            axs[-1 - j][i].set_xticks([])
            axs[-1 - j][i].set_ylim(bottom=0, top=1)
            axs[-1 - j][i].set_xlabel(chromosome[3:])
            axs[-1 - j][i].invert_yaxis()

            if i == 0:
                axs[-1 - j][i].set_ylabel(stage)

    matplotlib.pyplot.tight_layout()
    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
