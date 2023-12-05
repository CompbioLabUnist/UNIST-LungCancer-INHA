"""
aggregate_CNV_MutationSharedProportion_precancer.py: Aggregate CNV results w/ MSP as precancer vs. primary
"""
import argparse
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import tqdm
import step00

input_data = pandas.DataFrame()
watching = ""


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

    return numpy.average(a=a, weights=weights)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="CNV segment.tsv file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("size", help="SIZE file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--watching", help="Watching column name", type=str, required=True)
    parser.add_argument("--percentage", help="Percentage of patients to include", type=float, default=0.1)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0.0 < args.percentage < 0.5):
        raise ValueError("Percentage must be (0.0, 0.5)!!")

    watching = args.watching

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = list(clinical_data.index)
    print(len(patients), patients)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data = input_data.loc[(input_data["Patient"].isin(patients))]
    print(input_data)

    sample_list = sorted(set(input_data["Sample"]), key=step00.sorting)
    print(len(sample_list), sample_list)

    chromosome_list = list(filter(lambda x: x in set(input_data["chromosome"]), step00.chromosome_list))
    print(chromosome_list)

    size_data = pandas.read_csv(args.size, sep="\t", header=None, names=["chromosome", "length"]).set_index(keys="chromosome", verify_integrity=True)
    print(size_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    figures = list()
    for MSP in tqdm.tqdm(step00.sharing_columns):
        lower_bound, higher_bound = numpy.quantile(clinical_data[MSP], args.percentage), numpy.quantile(clinical_data[MSP], 1.0 - args.percentage)

        lower_precancer_list = list(clinical_data.loc[(clinical_data[MSP] <= lower_bound), f"{MSP}-sample"])
        higher_precancer_list = list(clinical_data.loc[(clinical_data[MSP] >= higher_bound), f"{MSP}-sample"])

        lower_primary_list = list(map(step00.get_paired_primary, lower_precancer_list))
        higher_primary_list = list(map(step00.get_paired_primary, higher_precancer_list))

        fig, axs = matplotlib.pyplot.subplots(ncols=len(chromosome_list), nrows=2, sharex="col", sharey="all", figsize=(64, 36), gridspec_kw={"width_ratios": list(map(lambda x: x / step00.big, size_data.loc[chromosome_list, "length"]))})

        for i, chromosome in enumerate(chromosome_list):
            chromosome_data = pandas.DataFrame(data=numpy.ones(shape=(len(lower_precancer_list) + len(lower_primary_list), size_data.loc[chromosome, "length"] // step00.big)), index=lower_precancer_list + lower_primary_list, dtype=float)

            with multiprocessing.Pool(args.cpus) as pool:
                for sample in lower_precancer_list + lower_primary_list:
                    chromosome_data.loc[sample, :] = pool.starmap(get_chromosome_data, [(sample, chromosome, step * step00.big, (step + 1) * step00.big) for step in list(chromosome_data.columns)])

            axs[0][i].errorbar(range(len(chromosome_data.columns)), numpy.mean(chromosome_data.loc[lower_primary_list, :], axis=0), yerr=numpy.std(chromosome_data.loc[lower_primary_list, :], axis=0), fmt="none", ecolor="gray", elinewidth=0.1)
            axs[0][i].errorbar(range(len(chromosome_data.columns)), numpy.mean(chromosome_data.loc[lower_precancer_list, :], axis=0), yerr=numpy.std(chromosome_data.loc[lower_precancer_list, :], axis=0) / 2, fmt="none", ecolor="tab:pink", elinewidth=0.1)

            axs[0][i].set_xticks([])
            axs[0][i].set_xlabel(chromosome[3:])
            axs[0][i].grid(True)

            if i == 0:
                axs[0][i].set_ylabel("Copy number")

        for i, chromosome in enumerate(chromosome_list):
            chromosome_data = pandas.DataFrame(data=numpy.ones(shape=(len(higher_precancer_list) + len(higher_primary_list), size_data.loc[chromosome, "length"] // step00.big)), index=higher_precancer_list + higher_primary_list, dtype=float)

            with multiprocessing.Pool(args.cpus) as pool:
                for sample in higher_precancer_list + higher_primary_list:
                    chromosome_data.loc[sample, :] = pool.starmap(get_chromosome_data, [(sample, chromosome, step * step00.big, (step + 1) * step00.big) for step in list(chromosome_data.columns)])

            axs[1][i].errorbar(range(len(chromosome_data.columns)), numpy.mean(chromosome_data.loc[higher_primary_list, :], axis=0), yerr=numpy.std(chromosome_data.loc[higher_primary_list, :], axis=0), fmt="none", ecolor="gray", elinewidth=0.1)
            axs[1][i].errorbar(range(len(chromosome_data.columns)), numpy.mean(chromosome_data.loc[higher_precancer_list, :], axis=0), yerr=numpy.std(chromosome_data.loc[higher_precancer_list, :], axis=0) / 2, fmt="none", ecolor="tab:pink", elinewidth=0.1)

            axs[1][i].set_xticks([])
            axs[1][i].set_xlabel(chromosome[3:])
            axs[1][i].grid(True)

            if i == 0:
                axs[1][i].set_ylabel("Copy number")

        matplotlib.pyplot.tight_layout()

        figures.append(f"{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

        print(MSP, "Lower", f"{lower_bound:.3f}", sorted(lower_precancer_list))
        print(MSP, "High", f"{higher_bound:.3f}", sorted(higher_precancer_list))

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure)
