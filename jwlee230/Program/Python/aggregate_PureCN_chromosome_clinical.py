"""
aggregate_PureCN_chromosome_clinical.py: Aggregate PureCN results as chromosome view with clincal information
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
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("size", help="SIZE file", type=str)
    parser.add_argument("cgc", help="CGC CSV files", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--compare", help="Comparison grouping (type, control, case)", type=str, nargs=3, default=["Recurrence", "NO", "YES"])
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

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="index", copy=False, ignore_index=True, verify_integrity=True)
    print(input_data)

    chromosome_list = list(filter(lambda x: x in set(input_data["chrom"]), step00.chromosome_full_list))

    control_sample_list = sorted(list(filter(lambda x: step00.get_patient(x) in control_patients, set(input_data["ID"]))), key=step00.sorting_by_type)
    control_primary_list = list(filter(lambda x: step00.get_long_sample_type(x) == "Primary", control_sample_list))
    control_precancer_list = list(filter(lambda x: step00.get_long_sample_type(x) != "Primary", control_sample_list))

    case_sample_list = sorted(list(filter(lambda x: step00.get_patient(x) in case_patients, set(input_data["ID"]))), key=step00.sorting_by_type)
    case_primary_list = list(filter(lambda x: step00.get_long_sample_type(x) == "Primary", case_sample_list))
    case_precancer_list = list(filter(lambda x: step00.get_long_sample_type(x) != "Primary", case_sample_list))

    size_data = pandas.read_csv(args.size, sep="\t", header=None, names=["chromosome", "length"]).set_index(keys="chromosome", verify_integrity=True)
    print(size_data)

    cgc_data = pandas.read_csv(args.cgc)
    cgc_data = cgc_data.loc[~(cgc_data["Genome Location"].str.contains(":-"))]
    with multiprocessing.Pool(args.cpus) as pool:
        cgc_data["Chromosome"] = pool.map(get_chromosome, cgc_data["Genome Location"])
        cgc_data["Start"] = pool.map(get_start, cgc_data["Genome Location"])
        cgc_data["End"] = pool.map(get_end, cgc_data["Genome Location"])
    print(cgc_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    for chromosome in tqdm.tqdm(chromosome_list):
        chromosome_data = pandas.DataFrame(data=numpy.ones(shape=(len(control_sample_list + case_sample_list), size_data.loc[chromosome, "length"] // step00.big)), index=control_sample_list + case_sample_list, dtype=float)

        for _, row in input_data.loc[(input_data["chrom"] == chromosome)].iterrows():
            chromosome_data.loc[row["ID"], row["loc.start"] // step00.big:row["loc.end"] // step00.big] = numpy.power(2, row[watching])

        fig, axs = matplotlib.pyplot.subplots(figsize=(64, 18 * 4), nrows=4, sharex=True)

        control_primary_proportion = list()
        control_precancer_proportion = list()
        case_primary_proportion = list()
        case_precancer_proportion = list()
        for j in range(chromosome_data.shape[1]):
            control_primary_proportion.append(len(list(filter(lambda x: chromosome_data.loc[x, j] >= (1 + args.threshold), control_primary_list))) / len(control_primary_list))
            control_precancer_proportion.append(len(list(filter(lambda x: chromosome_data.loc[x, j] >= (1 + args.threshold), control_precancer_list))) / len(control_precancer_list))
            case_primary_proportion.append(len(list(filter(lambda x: chromosome_data.loc[x, j] >= (1 + args.threshold), case_primary_list))) / len(case_primary_list))
            case_precancer_proportion.append(len(list(filter(lambda x: chromosome_data.loc[x, j] >= (1 + args.threshold), case_precancer_list))) / len(case_precancer_list))

        axs[0].plot(range(chromosome_data.shape[1]), control_primary_proportion, color="red", linestyle="--", label=args.compare[1])
        axs[0].plot(range(chromosome_data.shape[1]), case_primary_proportion, color="red", linestyle="-", label=args.compare[2])
        axs[0].set_ylim(bottom=0, top=1)
        axs[0].set_xlabel("{0} ({1:.1e} steps)".format(chromosome, step00.big))
        axs[0].set_ylabel("Primary")
        axs[0].legend(title=args.compare[0], loc="upper left")

        axs[1].plot(range(chromosome_data.shape[1]), control_precancer_proportion, color="lightsalmon", linestyle="--", label=args.compare[1])
        axs[1].plot(range(chromosome_data.shape[1]), case_precancer_proportion, color="lightsalmon", linestyle="-", label=args.compare[2])
        axs[1].set_ylim(bottom=0, top=1)
        axs[1].set_xlabel("{0} ({1:.1e} steps)".format(chromosome, step00.big))
        axs[1].set_ylabel("Precancer")
        axs[1].legend(title=args.compare[0], loc="upper left")

        control_primary_proportion = list()
        control_precancer_proportion = list()
        case_primary_proportion = list()
        case_precancer_proportion = list()
        for j in range(chromosome_data.shape[1]):
            control_primary_proportion.append(len(list(filter(lambda x: chromosome_data.loc[x, j] <= (1 - args.threshold), control_primary_list))) / len(control_primary_list))
            control_precancer_proportion.append(len(list(filter(lambda x: chromosome_data.loc[x, j] <= (1 - args.threshold), control_precancer_list))) / len(control_precancer_list))
            case_primary_proportion.append(len(list(filter(lambda x: chromosome_data.loc[x, j] <= (1 - args.threshold), case_primary_list))) / len(case_primary_list))
            case_precancer_proportion.append(len(list(filter(lambda x: chromosome_data.loc[x, j] <= (1 - args.threshold), case_precancer_list))) / len(case_precancer_list))

        axs[2].plot(range(chromosome_data.shape[1]), control_precancer_proportion, color="cyan", linestyle="--", label=args.compare[1])
        axs[2].plot(range(chromosome_data.shape[1]), case_precancer_proportion, color="cyan", linestyle="-", label=args.compare[2])
        axs[2].set_ylim(bottom=0, top=1)
        axs[2].invert_yaxis()
        axs[2].set_xlabel("{0} ({1:.1e} steps)".format(chromosome, step00.big))
        axs[2].set_ylabel("Precancer")
        axs[2].legend(title=args.compare[0], loc="lower left")

        axs[3].plot(range(chromosome_data.shape[1]), control_primary_proportion, color="navy", linestyle="--", label=args.compare[1])
        axs[3].plot(range(chromosome_data.shape[1]), case_primary_proportion, color="navy", linestyle="-", label=args.compare[2])
        axs[3].set_ylim(bottom=0, top=1)
        axs[3].invert_yaxis()
        axs[3].set_xlabel("{0} ({1:.1e} steps)".format(chromosome, step00.big))
        axs[3].set_ylabel("Priamry")
        axs[3].legend(title=args.compare[0], loc="lower left")

        texts = [[], [], [], []]
        for index, row in cgc_data.loc[(cgc_data["Chromosome"] == chromosome)].iterrows():
            x = (row["Start"] + row["End"]) / (2 * step00.big)

            for i in range(len(texts)):
                axs[i].axvline(x=x, color="k", linestyle="--")
                texts[i].append(axs[i].text(x, 0.5, row["Gene Symbol"], color="k", fontsize="small", horizontalalignment="left", bbox={"color": "white", "alpha": 0.8}))

        matplotlib.pyplot.tight_layout()
        for i in range(len(texts)):
            adjustText.adjust_text(texts[i], autoalign="y", ax=axs[i], lim=step00.big)
        figures.append("{0}.pdf".format(chromosome))
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
