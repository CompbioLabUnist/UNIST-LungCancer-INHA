"""
aggregate_signature_MSP_2.py: Aggregate mutational signature with MSP (relative)
"""
import argparse
import tarfile
import matplotlib
import matplotlib.colors
import matplotlib.patches
import matplotlib.pyplot
import numpy
import pandas
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("SBS", help="SBS input .TXT files", type=str)
    parser.add_argument("DBS", help="DBS input .TXT files", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--percentage", help="Percentage threshold", type=float, default=0.25)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".tsv"):
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
    print(len(patients), patients)

    SBS_data = pandas.read_csv(args.SBS, sep="\t", index_col=0)
    SBS_list = list(SBS_data.columns)
    SBS_data["Patient"] = list(map(step00.get_patient, list(SBS_data.index)))
    SBS_data = SBS_data.loc[(SBS_data["Patient"].isin(patients))]
    for index in tqdm.tqdm(list(SBS_data.index)):
        SBS_data.loc[index, SBS_list] = SBS_data.loc[index, SBS_list] / sum(SBS_data.loc[index, SBS_list])
    print(SBS_data)

    DBS_data = pandas.read_csv(args.DBS, sep="\t", index_col=0)
    DBS_list = list(DBS_data.columns)
    DBS_data["Patient"] = list(map(step00.get_patient, list(DBS_data.index)))
    DBS_data = DBS_data.loc[(DBS_data["Patient"].isin(patients))]
    for index in tqdm.tqdm(list(DBS_data.index)):
        DBS_data.loc[index, DBS_list] = DBS_data.loc[index, DBS_list] / sum(DBS_data.loc[index, DBS_list])
    print(DBS_data)

    SBS_color_dict = dict(zip(SBS_list, matplotlib.colors.TABLEAU_COLORS))
    DBS_color_dict = dict(zip(DBS_list, reversed(matplotlib.colors.TABLEAU_COLORS)))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    figures = list()
    for MSP in tqdm.tqdm(step00.sharing_columns[:2]):
        precancer_list = sorted(clinical_data[f"{MSP}-sample"], key=lambda x: clinical_data.loc[step00.get_patient(x), MSP])
        primary_list = list(map(step00.get_paired_primary, precancer_list))
        patient_list = list(map(step00.get_patient, precancer_list))

        MSP_Q1 = numpy.quantile(clinical_data[MSP], args.percentage)
        MSP_Q3 = numpy.quantile(clinical_data[MSP], 1.0 - args.percentage)

        mosaic = [["MSP-bar"], ["Survival"], ["SBS-Precancer"], ["DBS-Precancer"]]
        fig, axs = matplotlib.pyplot.subplot_mosaic(mosaic=mosaic, figsize=(18 * 3, 18 * 2), gridspec_kw={"height_ratios": [1, 1, 2, 2]}, layout="tight")

        MSP_Q_list = list(map(lambda x: "PSM-L" if (clinical_data.loc[x, MSP] <= MSP_Q1) else ("PSM-H" if (clinical_data.loc[x, MSP] >= MSP_Q3) else "None"), patient_list))
        MSP_L_list = list(map(lambda x: x[1], list(filter(lambda x: x[0] == "PSM-L", zip(MSP_Q_list, precancer_list)))))
        MSP_H_list = list(map(lambda x: x[1], list(filter(lambda x: x[0] == "PSM-H", zip(MSP_Q_list, precancer_list)))))

        axs["MSP-bar"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-L", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[step00.get_patient(x), MSP], MSP_L_list)), width=0.8, color="tab:blue", edgecolor=None, label="PSM-L")
        axs["MSP-bar"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-H", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[step00.get_patient(x), MSP], MSP_H_list)), width=0.8, color="tab:red", edgecolor=None, label="PSM-H")
        axs["MSP-bar"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "None", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[x[1], MSP], list(filter(lambda x: x[0] == "None", zip(MSP_Q_list, patient_list))))), width=0.8, color="tab:gray", edgecolor=None, label="Other")
        axs["MSP-bar"].set_xlabel("")
        axs["MSP-bar"].set_ylabel("PSM", fontsize="x-small")
        axs["MSP-bar"].set_xticks([])
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

        axs["Survival"].set_xticks([])
        axs["Survival"].set_yticks([0, 2500, 5000], ["0", "2500", "5000"], fontsize="xx-small", rotation="vertical", verticalalignment="center")
        axs["Survival"].set_ylabel("RFS (days)", fontsize="x-small")
        axs["Survival"].set_xlim(axs["MSP-bar"].get_xlim())
        axs["Survival"].legend(loc="upper left", fontsize="xx-small")
        axs["Survival"].grid(True)

        for i, signature in tqdm.tqdm(list(enumerate(SBS_list)), leave=False):
            axs["SBS-Precancer"].bar(x=range(len(precancer_list)), height=SBS_data.loc[precancer_list, :].iloc[:, i], bottom=SBS_data.loc[precancer_list, :].iloc[:, :i].sum(axis=1), color=SBS_color_dict[signature], edgecolor=None, linewidth=0, label=signature)
        axs["SBS-Precancer"].set_xticks(range(len(precancer_list)), ["" for _ in patient_list], fontsize="xx-small")
        axs["SBS-Precancer"].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"], fontsize="xx-small", rotation="vertical", verticalalignment="center")
        axs["SBS-Precancer"].set_ylabel("SBS in Precancer", fontsize="small")
        axs["SBS-Precancer"].set_ylim((0.0, 1.0))
        axs["SBS-Precancer"].grid(True)
        axs["SBS-Precancer"].legend(loc="upper right", fontsize="xx-small")

        for i, signature in tqdm.tqdm(list(enumerate(DBS_list)), leave=False):
            axs["DBS-Precancer"].bar(x=range(len(precancer_list)), height=DBS_data.loc[precancer_list, :].iloc[:, i], bottom=DBS_data.loc[precancer_list, :].iloc[:, :i].sum(axis=1), color=DBS_color_dict[signature], edgecolor=None, linewidth=0, label=signature)
        axs["DBS-Precancer"].set_xticks(range(len(precancer_list)), patient_list, fontsize="xx-small", rotation="vertical", horizontalalignment="center")
        axs["DBS-Precancer"].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"], fontsize="xx-small", rotation="vertical", verticalalignment="center")
        axs["DBS-Precancer"].set_ylabel("DBS in Precancer", fontsize="small")
        axs["DBS-Precancer"].set_ylim((0.0, 1.0))
        axs["DBS-Precancer"].grid(True)
        axs["DBS-Precancer"].legend(loc="upper right", fontsize="xx-small")

        figures.append(f"{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
