"""
aggregate_PureCN_chromosome.py: Aggregate PureCN results as chromosome view
"""
import argparse
import multiprocessing
import tarfile
import adjustText
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import tqdm
import step00

watching = "seg.mean"


def get_data(file_name: str) -> pandas.DataFrame:
    data = pandas.read_csv(file_name, sep="\t")
    data["chrom"] = list(map(lambda x: "chrX" if str(x) == "23" else ("chr" + str(x)), data["chrom"]))
    data["ID"] = step00.get_id(file_name)
    return data


def get_chromosome(location: str) -> str:
    return "chr" + location.split(":")[0]


def get_start(location: str) -> int:
    return int(location.replace("-", ":").split(":")[1])


def get_end(location: str) -> int:
    return int(location.replace("-", ":").split(":")[2])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PureCN output segments.TSV file(s)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("size", help="SIZE file", type=str)
    parser.add_argument("cgc", help="CGC CSV files", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)
    parser.add_argument("--threshold", help="Threshold for gain/loss", type=float, default=0.2)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.cgc.endswith(".csv"):
        raise ValueError("CGC must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
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

    sample_list = list(map(step00.get_id, args.input))

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="index", copy=False, ignore_index=True, verify_integrity=True)
    print(input_data)

    chromosome_list = list(filter(lambda x: x in set(input_data["chrom"]), step00.chromosome_full_list))
    print(chromosome_list)
    print(len(sample_list), sample_list)

    size_data = pandas.read_csv(args.size, sep="\t", header=None, names=["chromosome", "length"]).set_index(keys="chromosome", verify_integrity=True)
    print(size_data)

    cgc_data = pandas.read_csv(args.cgc)
    cgc_data = cgc_data.loc[~(cgc_data["Genome Location"].str.contains(":-"))]
    with multiprocessing.Pool(args.cpus) as pool:
        cgc_data["Chromosome"] = pool.map(get_chromosome, cgc_data["Genome Location"])
        cgc_data["Start"] = pool.map(get_start, cgc_data["Genome Location"])
        cgc_data["End"] = pool.map(get_end, cgc_data["Genome Location"])
    print(cgc_data)

    stage_set = set(map(step00.get_long_sample_type, args.input))
    stage_list = list(filter(lambda x: x in stage_set, step00.long_sample_type_list))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    for chromosome in tqdm.tqdm(chromosome_list):
        chromosome_data = pandas.DataFrame(data=numpy.ones(shape=(len(sample_list), size_data.loc[chromosome, "length"] // step00.big)), index=sample_list, dtype=float)

        for index, row in input_data.loc[(input_data["chrom"] == chromosome)].iterrows():
            chromosome_data.loc[row["ID"], row["loc.start"] // step00.big:row["loc.end"] // step00.big] = numpy.power(2, row[watching])

        fig, axs = matplotlib.pyplot.subplots(figsize=(64, 18 * 2), nrows=2, sharex=True)

        for j, stage in enumerate(stage_list):
            stage_sample_list = list(filter(lambda x: step00.get_long_sample_type(x) == stage, sample_list))
            proportion = [1 for _ in range(chromosome_data.shape[1])]
            for k in tqdm.tqdm(range(chromosome_data.shape[1])):
                proportion[k] = len(list(filter(lambda x: chromosome_data.loc[x, k] >= (1 + args.threshold), stage_sample_list))) / len(stage_sample_list)
            axs[0].plot(proportion, color=step00.stage_color_code[stage], linestyle=step00.stage_linestyle[stage], label=stage)

        axs[0].set_ylim(bottom=0, top=1)
        axs[0].set_xlabel("{0} ({1:.1e} steps)".format(chromosome, step00.big))
        axs[0].set_ylabel("Proportion")
        axs[0].legend(loc="upper left")

        for j, stage in enumerate(stage_list):
            stage_sample_list = list(filter(lambda x: step00.get_long_sample_type(x) == stage, sample_list))
            proportion = [1 for _ in range(chromosome_data.shape[1])]
            for k in tqdm.tqdm(range(chromosome_data.shape[1])):
                proportion[k] = len(list(filter(lambda x: chromosome_data.loc[x, k] <= (1 - args.threshold), stage_sample_list))) / len(stage_sample_list)
            axs[1].plot(proportion, color=step00.stage_color_code[stage], linestyle=step00.stage_linestyle[stage], label=stage)

        axs[1].set_ylim(bottom=0, top=1)
        axs[1].invert_yaxis()
        axs[1].set_xlabel("{0} ({1:.1e} steps)".format(chromosome, step00.big))
        axs[1].set_ylabel("Proportion")
        axs[1].legend(loc="lower left")

        upper_texts = list()
        lower_texts = list()
        for index, row in cgc_data.loc[(cgc_data["Chromosome"] == chromosome)].iterrows():
            x = (row["Start"] + row["End"]) / (2 * step00.big)

            axs[0].axvline(x=x, color="k", linestyle="--")
            axs[1].axvline(x=x, color="k", linestyle="--")

            upper_texts.append(axs[0].text(x, 0.5, row["Gene Symbol"], color="k", fontsize="small", horizontalalignment="center", bbox={"color": "white", "alpha": 0.8}))
            lower_texts.append(axs[1].text(x, 0.5, row["Gene Symbol"], color="k", fontsize="small", horizontalalignment="center", bbox={"color": "white", "alpha": 0.8}))

        matplotlib.pyplot.tight_layout()
        adjustText.adjust_text(upper_texts, autoalign="y", ax=axs[0], lim=step00.big)
        adjustText.adjust_text(lower_texts, autoalign="y", ax=axs[1], lim=step00.big)
        figures.append("{0}.pdf".format(chromosome))
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
