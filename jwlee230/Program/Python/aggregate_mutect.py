"""
aggregate_mutect.py: aggregate mutect MAF files
"""
import argparse
import collections
import multiprocessing
from comut import comut
import matplotlib
import numpy
import pandas
import step00


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("driver", help="MutEnricher Fisher enrichment output", type=str)
    parser.add_argument("census", help="Cancer gene census CSV file", type=str)
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.05)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.census.endswith(".csv"):
        raise ValueError("Census must end with .CSV!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.p < 1):
        raise ValueError("P-values must be (0, 1)")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    args.input = list(filter(lambda x: step00.get_patient(x.split("/")[-1].split(".")[0]) in patients, args.input))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    my_comut = comut.CoMut()
    my_comut.samples = sorted(list(map(lambda x: x.split("/")[-1].split(".")[0], args.input)), key=step00.sorting)

    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)

    mutect_data = mutect_data.loc[(mutect_data["Variant_Classification"].isin(["Nonsense_Mutation", "In_Frame_Del", "Frame_Shift_Ins", "Splice_Site", "In_Frame_Ins", "Frame_Shift_Del", "Missense_Mutation"]))]
    mutect_data["Tumor_Sample_Barcode"] = list(map(lambda x: x.split(".")[0], mutect_data["Tumor_Sample_Barcode"]))
    mutect_data["Variant_Classification"] = list(map(lambda x: {"Nonsense_Mutation": "Nonsense", "In_Frame_Del": "In frame indel", "In_Frame_Ins": "In frame indel", "Frame_Shift_Del": "Frameshift indel", "Missense_Mutation": "Missense", "Splice_Site": "Splice site", "Frame_Shift_Ins": "Frameshift indel"}[x], mutect_data["Variant_Classification"]))
    print(mutect_data)

    counter: collections.Counter = collections.Counter(mutect_data["Hugo_Symbol"])

    census_data = pandas.read_csv(args.census)
    print(census_data)

    driver_data = pandas.read_csv(args.driver, sep="\t")
    driver_data = driver_data.loc[(driver_data["Gene"].isin(mutect_data["Hugo_Symbol"])) & (driver_data["Gene"].isin(census_data["Gene Symbol"]))]
    for column in step00.MutEnricher_pval_columns:
        driver_data = driver_data.loc[(driver_data[column] < args.p)]
    driver_data.sort_values(by="Fisher_pval", ascending=False, ignore_index=True, inplace=True)
    driver_data["Count"] = list(map(lambda x: counter[x], driver_data["Gene"]))
    driver_data["-log10(P)"] = -1 * numpy.log10(driver_data["Fisher_pval"])
    print(driver_data)

    patient_data = pandas.DataFrame()
    patient_data["Tumor_Sample_Barcode"] = my_comut.samples
    patient_data["Patient"] = list(map(lambda x: hash(step00.get_patient(x)), my_comut.samples))
    patient_data["Collection_Type_category"] = "Type"
    patient_data["Collection_Type_value"] = list(map(step00.get_long_sample_type, my_comut.samples))
    patient_data["Mutation_Count"] = list(map(lambda x: mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == x)].shape[0], my_comut.samples))
    print(patient_data)

    my_comut.add_sample_indicators(patient_data[["Tumor_Sample_Barcode", "Patient"]].set_axis(labels=["sample", "group"], axis="columns"), name="Same patient")
    my_comut.add_categorical_data(patient_data[["Tumor_Sample_Barcode", "Collection_Type_category", "Collection_Type_value"]].set_axis(labels=["sample", "category", "value"], axis="columns"), name="Type", value_order=step00.long_sample_type_list)
    my_comut.add_categorical_data(mutect_data[["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification"]].set_axis(labels=["sample", "category", "value"], axis="columns"), name="Mutation type", category_order=driver_data["Gene"], mapping={"Missense": "green", "Nonsense": "deeppink", "In frame indel": {"facecolor": "blue", "hatch": "xxx"}, "Frameshift indel": "#FFD700", "Splice site": "darkviolet", "LOH": {"facecolor": "none", "edgecolor": "black", "linewidth": 3}, "Absent": {"facecolor": "grey", "alpha": 0.2}}, priority=["Frameshift indel"])
    my_comut.add_bar_data(patient_data[["Tumor_Sample_Barcode", "Mutation_Count"]].set_axis(labels=["sample", "group"], axis="columns"), name="Mutation count", ylabel="Counts", mapping={"group": "purple"})
    my_comut.add_side_bar_data(driver_data[["Gene", "-log10(P)"]].set_axis(labels=["category", "value"], axis="columns"), name="Mutation count", xlabel="-log10(P)", paired_name="Mutation type", position="left", mapping={"value": "tab:olive"})

    my_comut.plot_comut(x_padding=0.04, y_padding=0.04, tri_padding=0.03, figsize=(len(args.input), driver_data.shape[0] * 2))
    my_comut.add_unified_legend()
    my_comut.figure.savefig(args.output, bbox_inches="tight")
