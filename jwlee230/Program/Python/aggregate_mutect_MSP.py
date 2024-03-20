"""
aggregate_mutect_MSP.py: Aggregate mutect MAF files with MSP
"""
import argparse
import tarfile
import multiprocessing
from comut import comut
import matplotlib
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
    shared_data = shared_data.loc[(shared_data["Hugo_Symbol"].isin(gene_list)) & (shared_data["Variant_Classification"].isin(step00.nonsynonymous_mutations))]
    shared_data["Variant_Classification"] = "Shared"
    shared_data["Tumor_Sample_Barcode"] = shared_data["Precancer"]
    print(shared_data)

    mutect_data = pandas.concat([mutect_data, shared_data], join="inner", ignore_index=True)
    print(mutect_data)

    shared_data["Tumor_Sample_Barcode"] = shared_data["Primary"]
    mutect_data = pandas.concat([mutect_data, shared_data], join="inner", ignore_index=True)
    print(mutect_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    figures = list()
    for MSP in tqdm.tqdm(step00.sharing_columns):
        precancer_list = list(clinical_data[f"{MSP}-sample"])
        primary_list = list(map(step00.get_paired_primary, precancer_list))

        my_comut = comut.CoMut()
        my_comut.samples = sorted(precancer_list + primary_list, key=lambda x: clinical_data.loc[step00.get_patient(x), MSP])

        patient_data = pandas.DataFrame()
        patient_data["Tumor_Sample_Barcode"] = my_comut.samples
        patient_data["Patient"] = list(map(lambda x: hash(step00.get_patient(x)), my_comut.samples))
        patient_data["Collection_Type_category"] = "Type"
        patient_data["Collection_Type_value"] = list(map(lambda x: "Primary" if (step00.get_long_sample_type(x) == "Primary") else "Precancer", my_comut.samples))
        patient_data["MSP"] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), MSP], my_comut.samples))

        my_comut.add_sample_indicators(patient_data[["Tumor_Sample_Barcode", "Patient"]].set_axis(labels=["sample", "group"], axis="columns"), name="Same patient")
        my_comut.add_categorical_data(patient_data[["Tumor_Sample_Barcode", "Collection_Type_category", "Collection_Type_value"]].set_axis(labels=["sample", "category", "value"], axis="columns"), name="PRE/PRI", mapping={"Precancer": "tab:pink", "Primary": "gray"}, value_order=["Precancer", "Primary"])
        my_comut.add_categorical_data(mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"].isin(my_comut.samples)) & (mutect_data["Variant_Classification"].isin(set(step00.nonsynonymous_notations.values()) | {"Shared"})), ["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification"]].set_axis(labels=["sample", "category", "value"], axis="columns"), name="Mutation type", category_order=gene_list, priority=["Shared", "Frameshift indel"], mapping=step00.nonsynonymous_coloring | {"Shared": {"facecolor": "none", "edgecolor": "tab:red", "linewidth": 10}}, borders=["Shared"])
        my_comut.add_bar_data(patient_data[["Tumor_Sample_Barcode", "MSP"]].set_axis(labels=["sample", "MSP"], axis="columns"), name="MSP", ylabel="MSP", mapping={"MSP": "purple"})

        my_comut.plot_comut(x_padding=0.04, y_padding=0.04, tri_padding=0.03, figsize=(35 * 3, 35))
        my_comut.add_unified_legend()
        figures.append(f"{MSP}.pdf")
        my_comut.figure.savefig(figures[-1], bbox_inches="tight")

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
