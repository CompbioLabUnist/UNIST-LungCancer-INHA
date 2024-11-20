"""
aggregate_MSP_CNV_2.py: Aggregate CNV results with MSP as genome view
"""
import argparse
import collections
import multiprocessing
import tarfile
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import numpy
import pandas
import seaborn
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
        a.append(row[watching])
        weights.append((max(row["start"], start) - tmp_start + 1) / length)
        tmp_start = min(row["end"], end)

    a.append(center)
    weights.append((end - tmp_start + 1) / length)

    return numpy.average(a=a, weights=weights)


color_dict = {0: "navy", 1: "deepskyblue", 2: "lightgray", 3: "orange", 4: "chocolate", ">4": "darkred"}


def color_mapper(cnv):
    if (cnv < 0):
        raise ValueError("CNV value must be non-negative!!")

    cnv = round(cnv)

    if cnv > 4:
        return color_dict[">4"]
    else:
        return color_dict[cnv]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="CNV segment.tsv file", type=str)
    parser.add_argument("clinical", help="Clinical data TSV file", type=str)
    parser.add_argument("size", help="SIZE file", type=str)
    parser.add_argument("centromere", help="Centromere file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--watching", help="Watching column name", type=str, required=True)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)
    parser.add_argument("--percentage", help="Percentage threshold", type=float, default=0.25)
    parser.add_argument("--threshold", help="Threshold for gain/loss", type=float, default=0.2)

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
        raise ValueError("Percentage must be (0, 0.5)")
    elif not (0 < args.threshold < 1):
        raise ValueError("Threshold must be (0, 1)")

    watching = args.watching

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print("Patients:", len(patients))

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data = input_data.loc[(input_data["Patient"].isin(patients))]
    print(input_data)

    size_data = pandas.read_csv(args.size, sep="\t", header=None, names=["chromosome", "length"], index_col="chromosome")
    print(size_data)

    chromosome_list = list(filter(lambda x: x in set(input_data["chromosome"]), step00.chromosome_list))
    height_ratios = list(size_data.loc[chromosome_list, :].to_numpy() / 10 ** 7)

    centromere_data = pandas.read_csv(args.centromere, sep="\t")
    print(centromere_data)

    centromere_dict = dict()
    for chromosome in tqdm.tqdm(chromosome_list):
        centromere_dict[chromosome] = (min(centromere_data.loc[(centromere_data["chrom"] == chromosome), "chromStart"]) + max(centromere_data.loc[(centromere_data["chrom"] == chromosome), "chromEnd"])) / (2 * step00.big)
    print(centromere_dict)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    center = 2.0
    print("Min:", min(input_data[watching]))
    print("Max:", max(input_data[watching]))
    print(collections.Counter(input_data[watching]).most_common())

    figures = list()
    for MSP in tqdm.tqdm(step00.sharing_columns[:2]):
        precancer_list = sorted(clinical_data[f"{MSP}-sample"], key=lambda x: clinical_data.loc[step00.get_patient(x), MSP])
        primary_list = list(map(step00.get_paired_primary, precancer_list))
        patient_list = list(map(step00.get_patient, precancer_list))

        MSP_Q1 = numpy.quantile(clinical_data[MSP], args.percentage)
        MSP_Q3 = numpy.quantile(clinical_data[MSP], 1.0 - args.percentage)

        mosaic = [["MSP", "CNV-legend"], ["Survival", "CNV-legend"]] + [[chromosome, "."] for chromosome in chromosome_list]
        fig, axs = matplotlib.pyplot.subplot_mosaic(mosaic=mosaic, figsize=(18 * 3, 18 * 3), gridspec_kw={"height_ratios": [max(height_ratios), max(height_ratios)] + height_ratios, "width_ratios": [9.0, 1.0], "wspace": 0.0, "hspace": 0.0}, layout="tight")

        MSP_Q_list = list(map(lambda x: "PSM-L" if (clinical_data.loc[x, MSP] <= MSP_Q1) else ("PSM-H" if (clinical_data.loc[x, MSP] >= MSP_Q3) else "None"), patient_list))
        MSP_L_list = list(map(lambda x: x[1], list(filter(lambda x: x[0] == "PSM-L", zip(MSP_Q_list, precancer_list)))))
        MSP_H_list = list(map(lambda x: x[1], list(filter(lambda x: x[0] == "PSM-H", zip(MSP_Q_list, precancer_list)))))

        bar_list = list()
        bar_list.append(axs["MSP"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-L", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[step00.get_patient(x), MSP], MSP_L_list)), width=0.8, color="tab:blue", edgecolor=None, align="center", label="PSM-L"))
        bar_list.append(axs["MSP"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-H", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[step00.get_patient(x), MSP], MSP_H_list)), width=0.8, color="tab:red", edgecolor=None, align="center", label="PSM-H"))
        bar_list.append(axs["MSP"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "None", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[x[1], MSP], list(filter(lambda x: x[0] == "None", zip(MSP_Q_list, patient_list))))), width=0.8, color="tab:gray", edgecolor=None, align="center", label="Other"))

        axs["MSP"].set_xlabel("")
        axs["MSP"].set_ylabel("PSM", fontsize="x-small")
        axs["MSP"].set_xticks([])
        axs["MSP"].set_yticks([0.0, 0.25, 0.5], ["0", "0.25", "0.50"], fontsize="xx-small", rotation="vertical", verticalalignment="center")
        axs["MSP"].legend(loc="upper left", title="PSM", fontsize="xx-small")
        axs["MSP"].grid(True)

        survival_column = "Recurrence-Free Survival"
        axs["Survival"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-L", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[step00.get_patient(x), survival_column], MSP_L_list)), width=0.8, color="tab:blue", edgecolor=None, align="center")
        axs["Survival"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-H", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[step00.get_patient(x), survival_column], MSP_H_list)), width=0.8, color="tab:red", edgecolor=None, align="center")
        axs["Survival"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "None", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[x[1], survival_column], list(filter(lambda x: x[0] == "None", zip(MSP_Q_list, patient_list))))), width=0.8, color="tab:gray", edgecolor=None, align="center")

        recurrence_column = "Recurrence"
        axs["Survival"].scatter(list(filter(lambda x: clinical_data.loc[patient_list[x], recurrence_column] == "1", range(len(patient_list)))), list(map(lambda y: clinical_data.loc[patient_list[y], survival_column] + 500, list(filter(lambda x: clinical_data.loc[patient_list[x], recurrence_column] == "1", range(len(patient_list)))))), c="black", s=900, marker="*", edgecolor=None, label="Recurrence")

        axs["Survival"].set_xticks([])
        axs["Survival"].set_yticks([0, 2500, 5000], ["0", "2500", "5000"], fontsize="xx-small", rotation="vertical", verticalalignment="center")
        axs["Survival"].set_ylabel("RFS (days)", fontsize="x-small")
        axs["Survival"].set_xlim(axs["MSP"].get_xlim())
        axs["Survival"].set_ylim((0, 6000))
        axs["Survival"].legend(loc="upper left", fontsize="xx-small")
        axs["Survival"].grid(True)

        for chromosome in tqdm.tqdm(chromosome_list, leave=False):
            chromosome_data = pandas.DataFrame(data=numpy.ones(shape=(len(precancer_list), size_data.loc[chromosome, "length"] // step00.big)), index=precancer_list, dtype=float)
            precancer_data = pandas.DataFrame(index=patient_list, columns=range(size_data.loc[chromosome, "length"] // step00.big), dtype=str)

            with multiprocessing.Pool(args.cpus) as pool:
                for sample in precancer_list:
                    chromosome_data.loc[sample, :] = pool.starmap(get_chromosome_data, [(sample, chromosome, step * step00.big, (step + 1) * step00.big) for step in list(chromosome_data.columns)])

                for patient, precancer_sample in zip(patient_list, precancer_list):
                    precancer_data.loc[patient, :] = pool.map(color_mapper, chromosome_data.loc[precancer_sample, :])

            axs[chromosome].axhline(y=centromere_dict[chromosome], color="black", linewidth=5.0, linestyle="--")

            for i, patient in enumerate(patient_list):
                axs[chromosome].barh(y=list(chromosome_data.columns), width=[0.8 for _ in list(chromosome_data.columns)], left=i - 0.4, height=1.0, align="center", color=precancer_data.loc[patient, :], edgecolor=None, linewidth=0.0)
            axs[chromosome].set_ylabel(chromosome[3:], fontsize="xx-small")
            axs[chromosome].invert_yaxis()
            axs[chromosome].set_yticks([])
            axs[chromosome].set_xlim(axs["MSP"].get_xlim())
            axs[chromosome].grid(visible=True, axis="y")
            if chromosome == chromosome_list[-1]:
                axs[chromosome].set_xticks(range(len(precancer_list)), patient_list, fontsize="x-small", rotation="vertical")
            else:
                axs[chromosome].set_xticks(range(len(precancer_list)), ["" for _ in precancer_list])

        axs["CNV-legend"].legend(handles=[matplotlib.patches.Patch(color=value, label=key) for key, value in color_dict.items()], title="Copy number", loc="center", fontsize="x-small")
        axs["CNV-legend"].axis("off")

        figures.append(f"{MSP}-PRE.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    for MSP in tqdm.tqdm(step00.sharing_columns[:2]):
        precancer_list = sorted(clinical_data[f"{MSP}-sample"], key=lambda x: clinical_data.loc[step00.get_patient(x), MSP])
        primary_list = list(map(step00.get_paired_primary, precancer_list))
        patient_list = list(map(step00.get_patient, precancer_list))

        MSP_Q1 = numpy.quantile(clinical_data[MSP], args.percentage)
        MSP_Q3 = numpy.quantile(clinical_data[MSP], 1.0 - args.percentage)

        mosaic = [["MSP", "CNV-legend"], ["Survival", "CNV-legend"]] + [[chromosome, "."] for chromosome in chromosome_list]
        fig, axs = matplotlib.pyplot.subplot_mosaic(mosaic=mosaic, figsize=(18 * 3, 18 * 3), gridspec_kw={"height_ratios": [max(height_ratios), max(height_ratios)] + height_ratios, "width_ratios": [9.0, 1.0], "wspace": 0.0, "hspace": 0.0}, layout="tight")

        MSP_Q_list = list(map(lambda x: "PSM-L" if (clinical_data.loc[x, MSP] <= MSP_Q1) else ("PSM-H" if (clinical_data.loc[x, MSP] >= MSP_Q3) else "None"), patient_list))
        MSP_L_list = list(map(lambda x: x[1], list(filter(lambda x: x[0] == "PSM-L", zip(MSP_Q_list, precancer_list)))))
        MSP_H_list = list(map(lambda x: x[1], list(filter(lambda x: x[0] == "PSM-H", zip(MSP_Q_list, precancer_list)))))

        bar_list = list()
        bar_list.append(axs["MSP"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-L", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[step00.get_patient(x), MSP], MSP_L_list)), width=0.8, color="tab:blue", edgecolor=None, align="center", label="PSM-L"))
        bar_list.append(axs["MSP"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-H", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[step00.get_patient(x), MSP], MSP_H_list)), width=0.8, color="tab:red", edgecolor=None, align="center", label="PSM-H"))
        bar_list.append(axs["MSP"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "None", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[x[1], MSP], list(filter(lambda x: x[0] == "None", zip(MSP_Q_list, patient_list))))), width=0.8, color="tab:gray", edgecolor=None, align="center", label="Other"))

        axs["MSP"].set_xlabel("")
        axs["MSP"].set_ylabel("PSM", fontsize="x-small")
        axs["MSP"].set_xticks([])
        axs["MSP"].set_yticks([0.0, 0.25, 0.5], ["0", "0.25", "0.50"], fontsize="xx-small", rotation="vertical", verticalalignment="center")
        axs["MSP"].legend(loc="upper left", title="PSM", fontsize="xx-small")
        axs["MSP"].grid(True)

        survival_column = "Recurrence-Free Survival"
        axs["Survival"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-L", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[step00.get_patient(x), survival_column], MSP_L_list)), width=0.8, color="tab:blue", edgecolor=None, align="center")
        axs["Survival"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-H", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[step00.get_patient(x), survival_column], MSP_H_list)), width=0.8, color="tab:red", edgecolor=None, align="center")
        axs["Survival"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "None", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[x[1], survival_column], list(filter(lambda x: x[0] == "None", zip(MSP_Q_list, patient_list))))), width=0.8, color="tab:gray", edgecolor=None, align="center")

        recurrence_column = "Recurrence"
        axs["Survival"].scatter(list(filter(lambda x: clinical_data.loc[patient_list[x], recurrence_column] == "1", range(len(patient_list)))), list(map(lambda y: clinical_data.loc[patient_list[y], survival_column] + 500, list(filter(lambda x: clinical_data.loc[patient_list[x], recurrence_column] == "1", range(len(patient_list)))))), c="black", s=900, marker="*", edgecolor=None, label="Recurrence")

        axs["Survival"].set_xticks([])
        axs["Survival"].set_yticks([0, 2500, 5000], ["0", "2500", "5000"], fontsize="xx-small", rotation="vertical", verticalalignment="center")
        axs["Survival"].set_ylabel("RFS (days)", fontsize="x-small")
        axs["Survival"].set_xlim(axs["MSP"].get_xlim())
        axs["Survival"].set_ylim((0, 6000))
        axs["Survival"].legend(loc="upper left", fontsize="xx-small")
        axs["Survival"].grid(True)

        for chromosome in tqdm.tqdm(chromosome_list, leave=False):
            chromosome_data = pandas.DataFrame(data=numpy.ones(shape=(len(precancer_list), size_data.loc[chromosome, "length"] // step00.big)), index=precancer_list, dtype=float)
            primary_data = pandas.DataFrame(index=patient_list, columns=range(size_data.loc[chromosome, "length"] // step00.big), dtype=str)

            with multiprocessing.Pool(args.cpus) as pool:
                for sample in primary_list:
                    chromosome_data.loc[sample, :] = pool.starmap(get_chromosome_data, [(sample, chromosome, step * step00.big, (step + 1) * step00.big) for step in list(chromosome_data.columns)])

                for patient, primary_sample in zip(patient_list, primary_list):
                    primary_data.loc[patient, :] = pool.map(color_mapper, chromosome_data.loc[primary_sample, :])

            axs[chromosome].axhline(y=centromere_dict[chromosome], color="black", linewidth=5.0, linestyle="--")

            for i, patient in enumerate(patient_list):
                axs[chromosome].barh(y=list(chromosome_data.columns), width=[0.8 for _ in list(chromosome_data.columns)], left=i - 0.4, height=1.0, align="center", color=primary_data.loc[patient, :], edgecolor=None, linewidth=0.0)
            axs[chromosome].set_ylabel(chromosome[3:], fontsize="xx-small")
            axs[chromosome].invert_yaxis()
            axs[chromosome].set_yticks([])
            axs[chromosome].set_xlim(axs["MSP"].get_xlim())
            axs[chromosome].grid(visible=True, axis="y")
            if chromosome == chromosome_list[-1]:
                axs[chromosome].set_xticks(range(len(primary_list)), patient_list, fontsize="x-small", rotation="vertical")
            else:
                axs[chromosome].set_xticks(range(len(primary_list)), ["" for _ in primary_list])

        axs["CNV-legend"].legend(handles=[matplotlib.patches.Patch(color=value, label=key) for key, value in color_dict.items()], title="Copy number", loc="center", fontsize="x-small")
        axs["CNV-legend"].axis("off")

        figures.append(f"{MSP}-PRI.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
