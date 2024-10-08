"""
aggregate_mutect_cluster.py: aggregate mutect MAF files as cluster map
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

color_dict = {"Somatic": "tab:green", "Germline": "tab:pink", "Germline in PON": "tab:purple", "PON": "tab:cyan", "Other": "tab:gray"}


def check_quality(quality: str) -> str:
    if quality == "PASS":
        return "Somatic"
    elif ("normal_artifact" in quality) and ("panel_of_normals" in quality):
        return "Germline in PON"
    elif ("normal_artifact" in quality):
        return "Germline"
    elif ("panel_of_normals" in quality):
        return "PON"
    else:
        return "Other"


def read_vcf(filename: str) -> pandas.DataFrame:
    data = pandas.read_csv(filename, sep="\t", header=None, names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Normal", "Tumor"], comment="#")
    data["Sample"] = step00.get_id(filename)
    data["Quality"] = data["FILTER"].apply(check_quality)
    return data


def query(sample: str, quality: str) -> int:
    return len(mutect_data.loc[(mutect_data["Sample"] == sample) & (mutect_data["Quality"] == quality)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .VCF files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data TSV file w/ mutation shared proprotion", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--percentage", help="Percentage threshold", type=float, default=0.25)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".vcf"), args.input)):
        raise ValueError("INPUT must end with .VCF!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0.0 < args.percentage <= 0.5):
        raise ValueError("Percentage must be (0, 0.5]")

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(clinical_data)

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))
    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_vcf, args.input), ignore_index=True, copy=False)
    mutect_data = mutect_data.loc[(mutect_data["CHROM"].isin(step00.chromosome_list))]
    print(mutect_data)

    count_data = pandas.DataFrame(data=numpy.zeros((len(args.input), len(color_dict.keys()))), index=list(map(lambda x: step00.get_id(x), args.input)), columns=color_dict.keys(), dtype=int)
    with multiprocessing.Pool(args.cpus) as pool:
        for quality in tqdm.tqdm(list(color_dict.keys())):
            count_data[quality] = pool.starmap(query, [(sample, quality) for sample in list(map(lambda x: step00.get_id(x), args.input))])
    print(count_data)

    proportion_data = pandas.DataFrame(data=numpy.zeros((len(args.input), len(color_dict.keys()))), index=list(map(lambda x: step00.get_id(x), args.input)), columns=color_dict.keys(), dtype=float)
    for sample in tqdm.tqdm(list(map(lambda x: step00.get_id(x), args.input))):
        proportion_data.loc[sample, :] = count_data.loc[sample, :] / sum(count_data.loc[sample, :])
    print(proportion_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    figures = list()
    for MSP in tqdm.tqdm(step00.sharing_columns[:2]):
        precancer_list = sorted(clinical_data[f"{MSP}-sample"], key=lambda x: clinical_data.loc[step00.get_patient(x), MSP])
        primary_list = list(map(step00.get_paired_primary, precancer_list))
        patient_list = list(map(step00.get_patient, precancer_list))

        MSP_Q1 = numpy.quantile(clinical_data[MSP], args.percentage)
        MSP_Q3 = numpy.quantile(clinical_data[MSP], 1.0 - args.percentage)

        mosaic = [["MSP-bar"], ["Survival"], ["Count-PRE"], ["Count-PRI"], ["Proportion-PRE"], ["Proportion-PRI"]]
        fig, axs = matplotlib.pyplot.subplot_mosaic(mosaic=mosaic, figsize=(18 * 3, 18 * 4), gridspec_kw={"height_ratios": [1, 1, 2, 2, 2, 2]}, layout="tight")

        MSP_Q_list = list(map(lambda x: "PSM-L" if (clinical_data.loc[x, MSP] <= MSP_Q1) else ("PSM-H" if (clinical_data.loc[x, MSP] >= MSP_Q3) else "None"), patient_list))
        MSP_L_list = list(map(lambda x: x[1], list(filter(lambda x: x[0] == "PSM-L", zip(MSP_Q_list, precancer_list)))))
        MSP_H_list = list(map(lambda x: x[1], list(filter(lambda x: x[0] == "PSM-H", zip(MSP_Q_list, precancer_list)))))

        axs["MSP-bar"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-L", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[step00.get_patient(x), MSP], MSP_L_list)), width=0.8, color="tab:blue", edgecolor=None, label="PSM-L")
        axs["MSP-bar"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-H", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[step00.get_patient(x), MSP], MSP_H_list)), width=0.8, color="tab:red", edgecolor=None, label="PSM-H")
        axs["MSP-bar"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "None", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[x[1], MSP], list(filter(lambda x: x[0] == "None", zip(MSP_Q_list, patient_list))))), width=0.8, color="tab:gray", edgecolor=None, label="Other")
        axs["MSP-bar"].set_xlabel("")
        axs["MSP-bar"].set_ylabel("PSM", fontsize="x-small")
        axs["MSP-bar"].set_xticks(range(len(patient_list)), ["" for _ in patient_list])
        axs["MSP-bar"].set_yticks([0.0, 0.25, 0.5], ["0.00", "0.25", "0.50"], fontsize="xx-small", rotation="vertical", verticalalignment="center")
        axs["MSP-bar"].legend(loc="upper left", fontsize="xx-small")
        axs["MSP-bar"].grid(True)

        survival_column = "Recurrence-Free Survival"
        axs["Survival"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-L", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[step00.get_patient(x), survival_column], MSP_L_list)), width=0.8, color="tab:blue", edgecolor=None,)
        axs["Survival"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-H", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[step00.get_patient(x), survival_column], MSP_H_list)), width=0.8, color="tab:red", edgecolor=None)
        axs["Survival"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "None", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[x[1], survival_column], list(filter(lambda x: x[0] == "None", zip(MSP_Q_list, patient_list))))), width=0.8, color="tab:gray", edgecolor=None)

        recurrence_column = "Recurrence"
        recurrence_x_list = list(filter(lambda x: clinical_data.loc[patient_list[x], recurrence_column] == "1", range(len(patient_list))))
        axs["Survival"].scatter(recurrence_x_list, [1000 for _ in recurrence_x_list], c="black", s=1000, marker="*", edgecolor=None, label="Recurrence")

        axs["Survival"].set_xticks(range(len(patient_list)), ["" for _ in patient_list])
        axs["Survival"].set_yticks([0, 2500, 5000], ["0", "2500", "5000"], fontsize="xx-small", rotation="vertical", verticalalignment="center")
        axs["Survival"].set_ylabel("RFS (days)", fontsize="x-small")
        axs["Survival"].set_xlim(axs["MSP-bar"].get_xlim())
        axs["Survival"].legend(loc="upper left", fontsize="xx-small")
        axs["Survival"].grid(True)

        precancer_data = count_data.loc[precancer_list, :]
        for i, (quality, color) in tqdm.tqdm(list(enumerate(color_dict.items())), leave=False):
            axs["Count-PRE"].bar(x=range(len(precancer_list)), height=precancer_data.iloc[:, i], width=0.8, bottom=precancer_data.iloc[:, :i].sum(axis=1), color=color, edgecolor=None, align="center", label=quality)

        axs["Count-PRE"].set_xticks(range(len(patient_list)), ["" for _ in patient_list])
        axs["Count-PRE"].set_yticks([0, 50000, 100000, 150000, 200000, 250000], ["0", "50K", "100K", "150K", "200K", "250K"], fontsize="xx-small", rotation="vertical", verticalalignment="center")
        axs["Count-PRE"].set_ylabel("Mutation count in Precancer", fontsize="x-small")
        axs["Count-PRE"].set_xlim(axs["MSP-bar"].get_xlim())
        axs["Count-PRE"].legend(loc="upper left", fontsize="xx-small")
        axs["Count-PRE"].grid(True)

        primary_data = count_data.loc[primary_list, :]
        for i, (quality, color) in tqdm.tqdm(list(enumerate(color_dict.items())), leave=False):
            axs["Count-PRI"].bar(x=range(len(primary_list)), height=primary_data.iloc[:, i], width=0.8, bottom=primary_data.iloc[:, :i].sum(axis=1), color=color, edgecolor=None, align="center", label=quality)

        axs["Count-PRI"].set_xticks(range(len(patient_list)), ["" for _ in patient_list])
        axs["Count-PRI"].set_yticks([0, 50000, 100000, 150000, 200000, 250000], ["0", "50K", "100K", "150K", "200K", "250K"], fontsize="xx-small", rotation="vertical", verticalalignment="center")
        axs["Count-PRI"].set_ylabel("Mutation count in Primary", fontsize="x-small")
        axs["Count-PRI"].set_xlim(axs["MSP-bar"].get_xlim())
        axs["Count-PRI"].legend(loc="upper left", fontsize="xx-small")
        axs["Count-PRI"].grid(True)

        precancer_data = proportion_data.loc[precancer_list, :]
        for i, (quality, color) in tqdm.tqdm(list(enumerate(color_dict.items())), leave=False):
            axs["Proportion-PRE"].bar(x=range(len(precancer_list)), height=precancer_data.iloc[:, i], width=0.8, bottom=precancer_data.iloc[:, :i].sum(axis=1), color=color, edgecolor=None, align="center", label=quality)

        axs["Proportion-PRE"].set_xticks(range(len(patient_list)), ["" for _ in patient_list])
        axs["Proportion-PRE"].set_yticks([0, 0.5, 1], ["0.0", "0.5", "1.0"], fontsize="xx-small", rotation="vertical", verticalalignment="center")
        axs["Proportion-PRE"].set_xlim(axs["MSP-bar"].get_xlim())
        axs["Proportion-PRE"].set_ylabel("Mutation proportion in Precancer", fontsize="x-small")
        axs["Proportion-PRE"].set_ylim((0, 1))
        axs["Proportion-PRE"].legend(loc="upper left", fontsize="xx-small")
        axs["Proportion-PRE"].grid(True)

        primary_data = proportion_data.loc[primary_list, :]
        for i, (quality, color) in tqdm.tqdm(list(enumerate(color_dict.items())), leave=False):
            axs["Proportion-PRI"].bar(x=range(len(primary_list)), height=primary_data.iloc[:, i], width=0.8, bottom=primary_data.iloc[:, :i].sum(axis=1), color=color, edgecolor=None, align="center", label=quality)

        axs["Proportion-PRI"].set_xticks(range(len(patient_list)), patient_list, fontsize="xx-small", rotation="vertical", verticalalignment="center")
        axs["Proportion-PRI"].set_yticks([0, 0.5, 1], ["0.0", "0.5", "1.0"], fontsize="xx-small", rotation="vertical", verticalalignment="center")
        axs["Proportion-PRI"].set_xlim(axs["MSP-bar"].get_xlim())
        axs["Proportion-PRI"].set_ylabel("Mutation proportion in Primary", fontsize="x-small")
        axs["Proportion-PRI"].set_ylim((0, 1))
        axs["Proportion-PRI"].legend(loc="upper left", fontsize="xx-small")
        axs["Proportion-PRI"].grid(True)

        figures.append(f"{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
