"""
aggregate_mutect_MSP.py: Aggregate mutect MAF files with MSP
"""
import argparse
import tarfile
import multiprocessing
import matplotlib
import matplotlib.patches
import matplotlib.pyplot
import numpy
import pandas
import tqdm
import step00


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("driver", help="MutEnricher Fisher enrichment output", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("share", help="Mutation shared information TSV file", type=str)
    parser.add_argument("cgc", help="CGC CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.05)
    parser.add_argument("--percentage", help="Percentage threshold", type=float, default=0.25)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.share.endswith(".tsv"):
        raise ValueError("Share must end with .TSV!!")
    elif not args.cgc.endswith(".csv"):
        raise ValueError("CGC must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.p < 1):
        raise ValueError("P-values must be (0, 1)")
    elif not (0.0 < args.percentage < 0.5):
        raise ValueError("Percentage must be (0, 0.5)")

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

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))

    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
        mutect_data["Tumor_Sample_Barcode"] = pool.map(step00.get_id, mutect_data["Tumor_Sample_Barcode"])
    mutect_data["Variant_Classification"] = list(map(lambda x: step00.nonsynonymous_notations[x] if (x in step00.nonsynonymous_notations) else "Synonymous", mutect_data["Variant_Classification"]))
    mutect_data["Cancer_Stage"] = list(map(step00.get_long_sample_type, mutect_data["Tumor_Sample_Barcode"]))
    print(mutect_data)

    cgc_data = pandas.read_csv(args.cgc, index_col=0)
    print(cgc_data)

    driver_data = pandas.read_csv(args.driver, sep="\t")
    driver_data = driver_data.loc[(driver_data["Gene"].isin(mutect_data["Hugo_Symbol"])) & (driver_data["Gene"].isin(set(cgc_data.index)))]
    print(driver_data)

    for column in tqdm.tqdm(step00.MutEnricher_pval_columns):
        driver_data = driver_data.loc[(driver_data[column] < args.p)]

    driver_data.sort_values(by="Fisher_pval", ascending=False, inplace=True)
    print(driver_data)

    gene_list = list(driver_data["Gene"])[-40:]
    print("Gene:", len(gene_list))

    driver_data = driver_data.loc[(driver_data["Gene"].isin(gene_list))]
    mutect_data = mutect_data.loc[(mutect_data["Hugo_Symbol"].isin(gene_list))]
    print(mutect_data)

    shared_data = pandas.read_csv(args.share, sep="\t", index_col=0)
    shared_data = shared_data.loc[(shared_data["Hugo_Symbol"].isin(gene_list))]
    shared_data["Variant_Classification"] = "Shared"
    shared_data["Tumor_Sample_Barcode"] = shared_data["Precancer"]
    print(shared_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    figures = list()
    for MSP in tqdm.tqdm(step00.sharing_columns):
        precancer_list = sorted(clinical_data[f"{MSP}-sample"], key=lambda x: clinical_data.loc[step00.get_patient(x), MSP])
        primary_list = list(map(step00.get_paired_primary, precancer_list))
        patient_list = list(map(step00.get_patient, precancer_list))

        MSP_Q1 = numpy.quantile(clinical_data[MSP], args.percentage)
        MSP_Q3 = numpy.quantile(clinical_data[MSP], 1.0 - args.percentage)

        mosaic = [["MSP-bar", "Mutation-legend"], ["MSP-Q", "MSP-Q-legend"], ["TMB", "TMB-legend"], ["Age", "Age-legend"], ["Survival", "Survival-legend"], ["Stage", "Stage-legend"], ["Mutation", "Mutation-proportion"]]
        fig, axs = matplotlib.pyplot.subplot_mosaic(mosaic=mosaic, figsize=(18 * 3, 18 * 4), gridspec_kw={"height_ratios": [2, 1, 1, 1, 1, 1, 20], "width_ratios": [9, 1]}, layout="constrained")

        axs["MSP-bar"].bar(x=range(len(precancer_list)), height=clinical_data.loc[list(map(step00.get_patient, precancer_list)), MSP], width=0.8, color="violet", edgecolor=None)
        median_line = axs["MSP-bar"].axhline(y=numpy.median(clinical_data[MSP]), color="black", linestyle="--", linewidth=5)
        Q1_line = axs["MSP-bar"].axhline(y=MSP_Q1, color="black", linestyle=":", linewidth=5)
        Q3_line = axs["MSP-bar"].axhline(y=MSP_Q3, color="black", linestyle="-.", linewidth=5)
        axs["MSP-bar"].set_xlabel("")
        axs["MSP-bar"].set_ylabel("PSM")
        axs["MSP-bar"].set_xticks([])
        axs["MSP-bar"].legend(handles=[Q1_line, median_line, Q3_line], labels=["MSP-L", "Median", "MSP-H"], loc="upper left")

        MSP_Q_list = list(map(lambda x: "PSM-L" if (clinical_data.loc[x, MSP] <= MSP_Q1) else ("PSM-H" if (clinical_data.loc[x, MSP] >= MSP_Q3) else "None"), patient_list))
        axs["MSP-Q"].scatter(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-L", range(len(precancer_list)))), y=[0 for _ in list(filter(lambda x: MSP_Q_list[x] == "PSM-L", range(len(precancer_list))))], s=900, marker="o", c="tab:blue", edgecolor=None)
        axs["MSP-Q"].scatter(x=list(filter(lambda x: MSP_Q_list[x] == "None", range(len(precancer_list)))), y=[0 for _ in list(filter(lambda x: MSP_Q_list[x] == "None", range(len(precancer_list))))], s=900, marker="o", c="tab:gray", edgecolor=None)
        axs["MSP-Q"].scatter(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-H", range(len(precancer_list)))), y=[0 for _ in list(filter(lambda x: MSP_Q_list[x] == "PSM-H", range(len(precancer_list))))], s=900, marker="o", c="tab:red", edgecolor=None)
        axs["MSP-Q"].set_xticks([])
        axs["MSP-Q"].set_yticks([])
        axs["MSP-Q"].set_ylabel("PSM")

        axs["MSP-Q-legend"].legend(handles=[matplotlib.patches.Patch(color="tab:blue", label="PSM-L"), matplotlib.patches.Patch(color="tab:red", label="PSM-H"), matplotlib.patches.Patch(color="tab:gray", label="Other")], title="PSM", loc="center")
        axs["MSP-Q-legend"].axis("off")

        TMB_precancer = list(map(lambda x: mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == x)].shape[0], precancer_list))
        TMB_primary = list(map(lambda x: mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == x)].shape[0], primary_list))
        axs["TMB"].bar(x=range(len(precancer_list)), height=TMB_precancer, width=-0.4, align="edge", color="tab:pink")
        axs["TMB"].bar(x=range(len(precancer_list)), height=TMB_primary, width=0.4, align="edge", color="tab:gray")
        axs["TMB"].set_xticks([])
        axs["TMB"].set_ylabel("TMB")

        axs["TMB-legend"].legend(handles=[matplotlib.patches.Patch(color="tab:pink", label="Precancer"), matplotlib.patches.Patch(color="tab:gray", label="Primary")], title="PRE/PRI", loc="center")
        axs["TMB-legend"].axis("off")

        age_scatter = axs["Age"].scatter(x=range(len(precancer_list)), y=[0 for _ in precancer_list], s=900, c=clinical_data.loc[patient_list, "Age"], cmap="Wistia", edgecolor=None)
        axs["Age"].set_xticks([])
        axs["Age"].set_yticks([])
        axs["Age"].set_ylabel("Age")

        axs["Age-legend"].legend(*age_scatter.legend_elements(num=6), ncols=2, title="Age (years)", loc="center")
        axs["Age-legend"].axis("off")

        survival_scatter = axs["Survival"].scatter(x=range(len(precancer_list)), y=[0 for _ in precancer_list], s=900, c=clinical_data.loc[patient_list, "Overall Survival"], cmap="hot", edgecolor=None)
        axs["Survival"].set_xticks([])
        axs["Survival"].set_yticks([])
        axs["Survival"].set_ylabel("OS")

        axs["Survival-legend"].legend(*survival_scatter.legend_elements(num=6), ncols=2, title="OS (days)", loc="center")
        axs["Survival-legend"].axis("off")

        stage_colors = {0: "whitesmoke", 1: "lightgray", 2: "darkgray", 3: "dimgray", 4: "black"}
        axs["Stage"].scatter(x=range(len(precancer_list)), y=[0 for _ in precancer_list], s=900, marker="o", c=list(map(lambda x: stage_colors[clinical_data.loc[x, "Stage"]], patient_list)), edgecolor=None)
        axs["Stage"].scatter(x=range(len(precancer_list)), y=[1 for _ in precancer_list], s=900, marker="o", c=list(map(lambda x: stage_colors[clinical_data.loc[x, "pN"]], patient_list)), edgecolor=None)
        axs["Stage"].scatter(x=range(len(precancer_list)), y=[2 for _ in precancer_list], s=900, marker="o", c=list(map(lambda x: stage_colors[clinical_data.loc[x, "pT"]], patient_list)), edgecolor=None)
        axs["Stage"].set_xticks([])
        axs["Stage"].set_yticks([0, 1, 2], ["Stage", "pN", "pT"])
        axs["Stage"].grid(False)

        axs["Stage-legend"].legend(handles=[matplotlib.patches.Patch(color=value, label=key) for key, value in stage_colors.items()], title="Stage", loc="center", ncols=2)
        axs["Stage-legend"].axis("off")

        precancer_patch = None
        primary_patch = None
        both_patch = None
        shared_patch = None

        for i, precancer in tqdm.tqdm(enumerate(precancer_list), total=len(precancer_list), leave=False):
            for j, gene in enumerate(gene_list):
                is_precancer = not mutect_data.loc[(mutect_data["Hugo_Symbol"] == gene) & (mutect_data["Tumor_Sample_Barcode"] == precancer)].empty
                is_primary = not mutect_data.loc[(mutect_data["Hugo_Symbol"] == gene) & (mutect_data["Tumor_Sample_Barcode"] == step00.get_paired_primary(precancer))].empty
                is_shared = not shared_data.loc[(shared_data["Hugo_Symbol"] == gene) & (shared_data["Precancer"] == precancer) & (shared_data["Primary"] == step00.get_paired_primary(precancer))].empty

                if (not is_precancer) and (not is_primary):
                    continue

                color_dict = {False: "none", True: "dimgray"}
                tmp = axs["Mutation"].plot(i, j, fillstyle="left", marker="o", markersize=30, markerfacecolor="tab:red" if is_shared else color_dict[is_precancer], markerfacecoloralt="tab:red" if is_shared else color_dict[is_primary], markeredgecolor="tab:red" if is_shared else "lightgray", markeredgewidth=2.0)

                if (precancer_patch is None) and (is_precancer) and (not is_primary):
                    precancer_patch = tmp[0]
                elif (primary_patch is None) and (not is_precancer) and (is_primary):
                    primary_patch = tmp[0]
                elif (both_patch is None) and (is_precancer) and (is_primary):
                    both_patch = tmp[0]
                elif (shared_patch is None) and (is_shared):
                    shared_patch = tmp[0]

        axs["Mutation"].set_xticks(range(len(precancer_list)), patient_list, rotation="vertical")
        axs["Mutation"].set_yticks(range(len(gene_list)), gene_list)
        axs["Mutation"].grid(True)

        axs["Mutation-legend"].legend(handles=[precancer_patch, primary_patch, both_patch, shared_patch], labels=["Only PRE", "Only PRI", "Both, not shared", "Shared"], title="Mutation notation", loc="center")
        axs["Mutation-legend"].axis("off")

        axs["Mutation-proportion"].barh(y=range(len(gene_list)), width=list(map(lambda x: len(set(mutect_data.loc[(mutect_data["Hugo_Symbol"] == x) & (mutect_data["Tumor_Sample_Barcode"].isin(precancer_list)), "Tumor_Sample_Barcode"])) / len(precancer_list), gene_list)), height=0.4, align="edge", color="tab:pink", edgecolor=None)
        axs["Mutation-proportion"].barh(y=range(len(gene_list)), width=list(map(lambda x: len(set(mutect_data.loc[(mutect_data["Hugo_Symbol"] == x) & (mutect_data["Tumor_Sample_Barcode"].isin(primary_list)), "Tumor_Sample_Barcode"])) / len(primary_list), gene_list)), height=-0.4, align="edge", color="tab:gray", edgecolor=None)
        axs["Mutation-proportion"].set_xlabel("Proportion")
        axs["Mutation-proportion"].set_yticks([])

        figures.append(f"{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
