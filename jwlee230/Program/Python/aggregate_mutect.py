"""
aggregate_mutect.py: aggregate mutect MAF files
"""
import argparse
import itertools
import collections
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
    parser.add_argument("census", help="Cancer gene census CSV file", type=str)
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.05)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_including = parser.add_mutually_exclusive_group(required=True)
    group_including.add_argument("--include", help="Use included gene with CGC", action="store_true", default=False)
    group_including.add_argument("--exclude", help="Use included gene with CGC", action="store_true", default=False)

    group_sorting = parser.add_mutually_exclusive_group()
    group_sorting.add_argument("--patient", help="Sorting along patients", action="store_true", default=False)
    group_sorting.add_argument("--stage", help="Sorting along stages", action="store_true", default=False)
    group_sorting.add_argument("--gene", help="Sorting along genes", action="store_true", default=False)

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

    if not (args.patient or args.stage or args.gene):
        args.patient = True

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(sorted(patients))

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))

    if args.stage:
        args.input.sort(key=step00.sorting_by_type)
    elif args.patient:
        args.input.sort(key=step00.sorting)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
        mutect_data = mutect_data.loc[(mutect_data["Variant_Classification"].isin(step00.nonsynonymous_mutations))]
        mutect_data["Tumor_Sample_Barcode"] = pool.map(step00.get_id, mutect_data["Tumor_Sample_Barcode"])
    mutect_data["Variant_Classification"] = list(map(lambda x: step00.nonsynonymous_notations[x], mutect_data["Variant_Classification"]))
    print(mutect_data)

    counter: collections.Counter = collections.Counter(mutect_data.drop_duplicates(subset=["Hugo_Symbol", "Tumor_Sample_Barcode"])["Hugo_Symbol"])

    census_data = pandas.read_csv(args.census)
    print(census_data)

    driver_data = pandas.read_csv(args.driver, sep="\t")

    if args.include:
        driver_data = driver_data.loc[(driver_data["Gene"].isin(mutect_data["Hugo_Symbol"])) & (driver_data["Gene"].isin(census_data["Gene Symbol"]))]
    elif args.exclude:
        driver_data = driver_data.loc[(driver_data["Gene"].isin(mutect_data["Hugo_Symbol"])) & ~(driver_data["Gene"].isin(census_data["Gene Symbol"]))]
    else:
        raise Exception("Something went wrong!!")

    for column in tqdm.tqdm(step00.MutEnricher_pval_columns):
        driver_data = driver_data.loc[(driver_data[column] < args.p)]

    driver_data.sort_values(by="Fisher_pval", ascending=False, ignore_index=True, inplace=True)
    driver_data["Rate"] = list(map(lambda x: counter[x] / len(args.input), driver_data["Gene"]))
    driver_data["-log10(P)"] = -1 * numpy.log10(driver_data["Fisher_pval"])
    print(driver_data)

    if args.gene:
        sorting_data = pandas.DataFrame(data=numpy.zeros((len(args.input), driver_data.shape[0])), index=args.input, columns=reversed(driver_data["Gene"]), dtype=bool)

        for patient, gene in tqdm.tqdm(list(itertools.product(list(sorting_data.index), list(sorting_data.columns)))):
            if not mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == step00.get_id(patient)) & (mutect_data["Hugo_Symbol"] == gene)].empty:
                sorting_data.loc[patient, gene] = True

        sorting_data.sort_values(by=list(sorting_data.columns), axis="index", ascending=False, kind="mergesort", inplace=True)
        args.input = list(sorting_data.index)

    my_comut = comut.CoMut()
    my_comut.samples = list(map(lambda x: x.split("/")[-1].split(".")[0], args.input))

    patient_data = pandas.DataFrame()
    patient_data["Tumor_Sample_Barcode"] = my_comut.samples
    patient_data["Patient"] = list(map(lambda x: hash(step00.get_patient(x)), my_comut.samples))
    patient_data["Collection_Type_category"] = "Type"
    patient_data["Collection_Type_value"] = list(map(step00.get_long_sample_type, my_comut.samples))
    patient_data["non-synonymous"] = list(map(lambda x: mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == x)].shape[0], my_comut.samples))
    print(patient_data)

    if args.patient:
        my_comut.add_sample_indicators(patient_data[["Tumor_Sample_Barcode", "Patient"]].set_axis(labels=["sample", "group"], axis="columns"), name="Same patient")

    my_comut.add_categorical_data(patient_data[["Tumor_Sample_Barcode", "Collection_Type_category", "Collection_Type_value"]].set_axis(labels=["sample", "category", "value"], axis="columns"), name="Stage" ,mapping=step00.stage_color_code)
    my_comut.add_categorical_data(mutect_data[["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification"]].set_axis(labels=["sample", "category", "value"], axis="columns"), name="Mutation type", category_order=driver_data["Gene"], priority=["Frameshift indel"], mapping=step00.nonsynonymous_coloring)
    my_comut.add_bar_data(patient_data[["Tumor_Sample_Barcode", "non-synonymous"]].set_axis(labels=["sample", "Counts"], axis="columns"), name="Mutation count", ylabel="TMB", mapping={"Counts": "purple"})
    # my_comut.add_side_bar_data(driver_data[["Gene", "-log10(P)"]].set_axis(labels=["category", "value"], axis="columns"), name="P-value", xlabel="-log10(P)", paired_name="Mutation type", position="left", mapping={"value": "olive"})
    my_comut.add_side_bar_data(driver_data[["Gene", "Rate"]].set_axis(labels=["category", "value"], axis="columns"), name="Mutation rate", xlabel="Present Rate", paired_name="Mutation type", position="left", mapping={"value": "teal"})

    my_comut.plot_comut(x_padding=0.04, y_padding=0.04, tri_padding=0.03, figsize=(len(args.input) + 10, driver_data.shape[0] * 2))
    my_comut.add_unified_legend()
    my_comut.figure.savefig(args.output, bbox_inches="tight")
