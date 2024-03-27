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

    cgc_data = pandas.read_csv(args.cgc, index_col=0).dropna(axis="index", subset=["Tumour Types(Somatic)"])
    cgc_data = cgc_data.loc[(cgc_data["Tumour Types(Somatic)"].str.contains("lung")) | (cgc_data["Tumour Types(Somatic)"].str.contains("NSCLC")) | (cgc_data["Tumour Types(Somatic)"].str.contains("ALL"))]
    print(cgc_data)

    driver_data = pandas.read_csv(args.driver, sep="\t")
    driver_data = driver_data.loc[(driver_data["Gene"].isin(set(mutect_data["Hugo_Symbol"]))) & (driver_data["Gene"].isin(set(cgc_data.index)))]
    print(driver_data)

    for column in tqdm.tqdm(step00.MutEnricher_pval_columns):
        driver_data = driver_data.loc[(driver_data[column] < args.p)]

    driver_data.sort_values(by="Fisher_pval", ascending=False, inplace=True)
    print(driver_data)

    gene_list = list(driver_data["Gene"])
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

        mosaic = [["MSP-bar", "Mutation-legend"], ["TMB", "TMB-legend"], ["Survival", "Mutation-proportion-legend"], ["Mutation", "Mutation-proportion"]]
        fig, axs = matplotlib.pyplot.subplot_mosaic(mosaic=mosaic, figsize=(18 * 3, 18 * 4), gridspec_kw={"height_ratios": [1, 1, 1, 20], "width_ratios": [9, 1]}, layout="tight")

        MSP_Q_list = list(map(lambda x: "PSM-L" if (clinical_data.loc[x, MSP] <= MSP_Q1) else ("PSM-H" if (clinical_data.loc[x, MSP] >= MSP_Q3) else "None"), patient_list))
        MSP_L_list = list(map(lambda x: x[1], list(filter(lambda x: x[0] == "PSM-L", zip(MSP_Q_list, precancer_list)))))
        MSP_H_list = list(map(lambda x: x[1], list(filter(lambda x: x[0] == "PSM-H", zip(MSP_Q_list, precancer_list)))))

        axs["MSP-bar"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-L", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[step00.get_patient(x), MSP], MSP_L_list)), width=0.8, color="tab:blue", edgecolor=None, label="PSM-L")
        axs["MSP-bar"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-H", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[step00.get_patient(x), MSP], MSP_H_list)), width=0.8, color="tab:red", edgecolor=None, label="PSM-H")
        axs["MSP-bar"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "None", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[x[1], MSP], list(filter(lambda x: x[0] == "None", zip(MSP_Q_list, patient_list))))), width=0.8, color="tab:gray", edgecolor=None, label="Other")
        axs["MSP-bar"].set_xlabel("")
        axs["MSP-bar"].set_ylabel("PSM")
        axs["MSP-bar"].set_xticks([])
        axs["MSP-bar"].set_yticks([0.0, 0.25, 0.5], ["0.0", "0.25", "0.5"], fontsize="xx-small", rotation="vertical")
        axs["MSP-bar"].legend(loc="upper left")
        axs["MSP-bar"].grid(True)

        TMB_precancer = list(map(lambda x: mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == x)].shape[0] / step00.WES_length * step00.big, precancer_list))
        TMB_primary = list(map(lambda x: mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == x)].shape[0] / step00.WES_length * step00.big, primary_list))
        axs["TMB"].bar(x=range(len(precancer_list)), height=TMB_precancer, width=-0.4, align="edge", color="tab:green", label="Precancer")
        axs["TMB"].bar(x=range(len(precancer_list)), height=TMB_primary, width=0.4, align="edge", color="tab:purple", label="Primary")
        axs["TMB"].set_xticks([])
        axs["TMB"].set_yticks([0.0, 0.5, 1.0, 1.5], ["0.0", "0.5", "1.0", "1.5"], fontsize="xx-small", rotation="vertical")
        axs["TMB"].set_ylabel("TMB (#/Mb)")
        axs["TMB"].set_xlim(axs["MSP-bar"].get_xlim())
        axs["TMB"].grid(True)

        axs["TMB-legend"].legend(handles=[matplotlib.patches.Patch(color=value, label=key) for key, value in step00.precancer_color_code.items()], title="PRE/PRI", loc="center")
        axs["TMB-legend"].axis("off")

        survival_column = "Overall Survival"
        axs["Survival"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-L", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[step00.get_patient(x), survival_column], MSP_L_list)), width=0.8, color="tab:blue", edgecolor=None, label="PSM-L")
        axs["Survival"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "PSM-H", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[step00.get_patient(x), survival_column], MSP_H_list)), width=0.8, color="tab:red", edgecolor=None, label="PSM-H")
        axs["Survival"].bar(x=list(filter(lambda x: MSP_Q_list[x] == "None", range(len(precancer_list)))), height=list(map(lambda x: clinical_data.loc[x[1], survival_column], list(filter(lambda x: x[0] == "None", zip(MSP_Q_list, patient_list))))), width=0.8, color="tab:gray", edgecolor=None, label="Other")
        axs["Survival"].set_xticks([])
        axs["Survival"].set_yticks([0, 2000, 4000], ["0", "2000", "4000"], fontsize="xx-small", rotation="vertical")
        axs["Survival"].set_ylabel("OS (days)")
        axs["Survival"].set_xlim(axs["MSP-bar"].get_xlim())
        axs["Survival"].legend(loc="upper left")
        axs["Survival"].grid(True)

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
                if is_shared:
                    tmp = axs["Mutation"].plot(i, j, fillstyle="left", marker="o", markersize=30, markerfacecolor="tab:olive", markerfacecoloralt="tab:olive", markeredgecolor="tab:olive" if is_shared else "lightgray", markeredgewidth=2.0)
                else:
                    tmp = axs["Mutation"].plot(i, j, fillstyle="left", marker="o", markersize=30, markerfacecolor=color_dict[is_precancer], markerfacecoloralt=color_dict[is_primary], markeredgecolor="lightgray", markeredgewidth=2.0)

                if (precancer_patch is None) and (is_precancer) and (not is_primary):
                    precancer_patch = tmp[0]
                elif (primary_patch is None) and (not is_precancer) and (is_primary):
                    primary_patch = tmp[0]
                elif (both_patch is None) and (is_precancer) and (is_primary):
                    both_patch = tmp[0]
                elif (shared_patch is None) and (is_shared):
                    shared_patch = tmp[0]

        axs["Mutation"].set_xticks(range(len(precancer_list)), patient_list, rotation="vertical")
        axs["Mutation"].set_yticks(range(len(gene_list)), gene_list, fontsize="small")
        axs["Mutation"].set_xlim(axs["MSP-bar"].get_xlim())
        axs["Mutation"].grid(True)

        axs["Mutation-legend"].legend(handles=[precancer_patch, primary_patch, both_patch, shared_patch], labels=["Only Precancer", "Only Primary", "Both, not shared", "Shared"], title="Mutation notation", loc="center")
        axs["Mutation-legend"].axis("off")

        bar_list = list()
        bar_list.append(axs["Mutation-proportion"].barh(y=range(len(gene_list)), width=list(map(lambda x: len(set(mutect_data.loc[(mutect_data["Hugo_Symbol"] == x) & (mutect_data["Tumor_Sample_Barcode"].isin(MSP_L_list)), "Tumor_Sample_Barcode"])) / len(MSP_L_list), gene_list)), height=0.4, align="edge", color="tab:cyan", edgecolor=None, label="PSM-L"))
        bar_list.append(axs["Mutation-proportion"].barh(y=range(len(gene_list)), width=list(map(lambda x: len(set(shared_data.loc[(shared_data["Hugo_Symbol"] == x) & (shared_data["Precancer"].isin(MSP_L_list)), "Tumor_Sample_Barcode"])) / len(MSP_L_list), gene_list)), height=0.4, align="edge", color="tab:blue", edgecolor=None, label="PSM-L (Shared)"))
        bar_list.append(axs["Mutation-proportion"].barh(y=range(len(gene_list)), width=list(map(lambda x: len(set(mutect_data.loc[(mutect_data["Hugo_Symbol"] == x) & (mutect_data["Tumor_Sample_Barcode"].isin(MSP_H_list)), "Tumor_Sample_Barcode"])) / len(MSP_H_list), gene_list)), height=-0.4, align="edge", color="tab:pink", edgecolor=None, label="PSM-H"))
        bar_list.append(axs["Mutation-proportion"].barh(y=range(len(gene_list)), width=list(map(lambda x: len(set(shared_data.loc[(shared_data["Hugo_Symbol"] == x) & (shared_data["Precancer"].isin(MSP_H_list)), "Tumor_Sample_Barcode"])) / len(MSP_L_list), gene_list)), height=-0.4, align="edge", color="tab:red", edgecolor=None, label="PSM-H (Shared)"))

        axs["Mutation-proportion"].set_xlabel("Proportion")
        axs["Mutation-proportion"].set_yticks([])
        axs["Mutation-proportion"].set_ylim(axs["Mutation"].get_ylim())
        axs["Mutation-proportion"].grid(True)

        axs["Mutation-proportion-legend"].legend(handles=bar_list, title="Mutations in PRE", loc="center")
        axs["Mutation-proportion-legend"].axis("off")

        figures.append(f"{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
