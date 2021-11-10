"""
aggregate_mutect_2.py: aggregate mutect MAF files with separation of Recur vs. Non-recur
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
    parser.add_argument("--compare", help="Comparison grouping (type, control, case)", type=str, nargs=3, default=["Recurrence", "NO", "YES"])
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.05)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

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
        case_patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC") & (clinical_data[args.compare[0]] == args.compare[2])].index)
        control_patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC") & (clinical_data[args.compare[0]] == args.compare[1])].index)
    elif args.ADC:
        case_patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC") & (clinical_data[args.compare[0]] == args.compare[2])].index)
        control_patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC") & (clinical_data[args.compare[0]] == args.compare[1])].index)
    else:
        raise Exception("Something went wrong!!")
    print(sorted(control_patients))
    print(sorted(case_patients))

    args.input = sorted(filter(lambda x: step00.get_patient(x.split("/")[-1].split(".")[0]) in control_patients, args.input), key=step00.sorting_by_type) + sorted(filter(lambda x: step00.get_patient(x.split("/")[-1].split(".")[0]) in case_patients, args.input), key=step00.sorting_by_type)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)

    mutect_data = mutect_data.loc[(mutect_data["Variant_Classification"].isin(step00.nonsynonymous_mutations))]
    mutect_data["Tumor_Sample_Barcode"] = list(map(lambda x: x.split(".")[0], mutect_data["Tumor_Sample_Barcode"]))
    mutect_data["Variant_Classification"] = list(map(lambda x: step00.nonsynonymous_notations[x], mutect_data["Variant_Classification"]))
    print(mutect_data)

    counter: collections.Counter = collections.Counter(mutect_data.drop_duplicates(subset=["Hugo_Symbol", "Tumor_Sample_Barcode"])["Hugo_Symbol"])

    census_data = pandas.read_csv(args.census)
    print(census_data)

    driver_data = pandas.read_csv(args.driver, sep="\t")
    driver_data = driver_data.loc[(driver_data["Gene"].isin(mutect_data["Hugo_Symbol"])) & (driver_data["Gene"].isin(census_data["Gene Symbol"]))]
    for column in step00.MutEnricher_pval_columns:
        driver_data = driver_data.loc[(driver_data[column] < args.p)]
    driver_data.sort_values(by="Fisher_pval", ascending=False, ignore_index=True, inplace=True)
    driver_data["Rate"] = list(map(lambda x: counter[x] / len(args.input), driver_data["Gene"]))
    driver_data["-log10(P)"] = -1 * numpy.log10(driver_data["Fisher_pval"])
    print(driver_data)

    my_comut = comut.CoMut()
    my_comut.samples = list(map(lambda x: x.split("/")[-1].split(".")[0], args.input))

    patient_data = pandas.DataFrame()
    patient_data["Tumor_Sample_Barcode"] = my_comut.samples
    patient_data["Patient"] = list(map(lambda x: hash(step00.get_patient(x)), my_comut.samples))
    patient_data["Collection_Type_category"] = "Type"
    patient_data["Collection_Type_value"] = list(map(step00.get_long_sample_type, my_comut.samples))
    patient_data[args.compare[0]] = list(map(lambda x: args.compare[2] if (step00.get_patient(x) in case_patients) else args.compare[1], my_comut.samples))
    patient_data["non-synonymous"] = list(map(lambda x: mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == x)].shape[0], my_comut.samples))
    print(patient_data)

    my_comut.add_categorical_data(patient_data[["Tumor_Sample_Barcode", "Collection_Type_category", "Collection_Type_value"]].set_axis(labels=["sample", "category", "value"], axis="columns"), name="Stage")
    my_comut.add_categorical_data(patient_data[["Tumor_Sample_Barcode", "Collection_Type_category", "Recur?"]].set_axis(labels=["sample", "category", "value"], axis="columns"), name="Recur?", mapping={"Recur": "tab:red", "Non-Recur": "tab:green"})
    my_comut.add_categorical_data(mutect_data[["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification"]].set_axis(labels=["sample", "category", "value"], axis="columns"), name="Mutation type", category_order=driver_data["Gene"], priority=["Frameshift indel"], mapping=step00.nonsynonymous_coloring)
    my_comut.add_bar_data(patient_data[["Tumor_Sample_Barcode", "non-synonymous"]].set_axis(labels=["sample", "Counts"], axis="columns"), name="Mutation count", ylabel="Mutations", mapping={"Counts": "purple"})
    my_comut.add_side_bar_data(driver_data[["Gene", "-log10(P)"]].set_axis(labels=["category", "value"], axis="columns"), name="P-value", xlabel="-log10(P)", paired_name="Mutation type", position="left", mapping={"value": "olive"})
    my_comut.add_side_bar_data(driver_data[["Gene", "Rate"]].set_axis(labels=["category", "value"], axis="columns"), name="Mutation rate", xlabel="Present Rate", paired_name="Mutation type", position="right", mapping={"value": "teal"})

    my_comut.plot_comut(x_padding=0.04, y_padding=0.04, tri_padding=0.03, figsize=(len(args.input), driver_data.shape[0] * 3))
    my_comut.add_unified_legend()
    my_comut.figure.savefig(args.output, bbox_inches="tight")
