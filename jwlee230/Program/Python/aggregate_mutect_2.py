"""
aggregate_mutect_2.py: aggregate mutect MAF files with separation of clinical data
"""
import argparse
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
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output file", type=str)
    parser.add_argument("--compare", help="Comparison grouping (type, control, case)", type=str, nargs="+", default=["Recurrence", "NO", "YES"])
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
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC") & (clinical_data[args.compare[0]].isin(args.compare[1:]))].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC") & (clinical_data[args.compare[0]].isin(args.compare[1:]))].index)
    else:
        raise Exception("Something went wrong!!")

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))
    tmp = list()
    for c in tqdm.tqdm(args.compare[1:]):
        tmp += sorted(filter(lambda x: (step00.get_patient(x) in set(clinical_data.loc[clinical_data[args.compare[0]] == c].index)), args.input), key=step00.sorting_by_type)
    args.input = tmp

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
        mutect_data["Tumor_Sample_Barcode"] = pool.map(step00.get_id, mutect_data["Tumor_Sample_Barcode"])
    mutect_data["Variant_Classification"] = list(map(lambda x: step00.nonsynonymous_notations[x] if (x in step00.nonsynonymous_notations) else "Synonymous", mutect_data["Variant_Classification"]))
    mutect_data["Cancer_Stage"] = list(map(step00.get_long_sample_type, mutect_data["Tumor_Sample_Barcode"]))
    print(mutect_data)

    census_data = pandas.read_csv(args.census)
    print(census_data)

    driver_data = pandas.read_csv(args.driver, sep="\t")
    driver_data = driver_data.loc[(driver_data["Gene"].isin(mutect_data["Hugo_Symbol"])) & (driver_data["Gene"].isin(census_data["Gene Symbol"]))]

    for column in tqdm.tqdm(step00.MutEnricher_pval_columns):
        driver_data = driver_data.loc[(driver_data[column] < args.p)]

    driver_data.sort_values(by="Fisher_pval", ascending=False, ignore_index=True, inplace=True)

    stage_list = list(filter(lambda x: x in set(mutect_data["Cancer_Stage"]), reversed(step00.long_sample_type_list)))
    print(stage_list)
    for stage in tqdm.tqdm(stage_list):
        counter: collections.Counter = collections.Counter(mutect_data.loc[mutect_data["Cancer_Stage"] == stage].drop_duplicates(subset=["Hugo_Symbol", "Tumor_Sample_Barcode"])["Hugo_Symbol"])
        driver_data[stage] = list(map(lambda x: counter[x] / len(args.input), driver_data["Gene"]))
    driver_data["-log10(P)"] = -1 * numpy.log10(driver_data["Fisher_pval"])
    print(driver_data)

    driver_data["-log10(P)"] = -1 * numpy.log10(driver_data["Fisher_pval"])
    print(driver_data)

    my_comut = comut.CoMut()
    my_comut.samples = list(map(step00.get_id, args.input))

    patient_data = pandas.DataFrame()
    patient_data["Tumor_Sample_Barcode"] = my_comut.samples
    patient_data["Patient"] = list(map(lambda x: hash(step00.get_patient(x)), my_comut.samples))
    patient_data["Collection_Type_category"] = "Type"
    patient_data["Collection_Type_value"] = list(map(step00.get_long_sample_type, my_comut.samples))
    patient_data[args.compare[0]] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), args.compare[0]], my_comut.samples))
    patient_data["TMB"] = list(map(lambda x: mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == x)].shape[0] / step00.WES_length * 10 ** 6, my_comut.samples))
    print(patient_data)

    my_comut.add_categorical_data(patient_data[["Tumor_Sample_Barcode", "Collection_Type_category", "Collection_Type_value"]].set_axis(labels=["sample", "category", "value"], axis="columns"), name="Stage", mapping=step00.stage_color_code, value_order=reversed(stage_list))
    my_comut.add_categorical_data(patient_data[["Tumor_Sample_Barcode", "Collection_Type_category", args.compare[0]]].set_axis(labels=["sample", "category", "value"], axis="columns"), name=args.compare[0])
    my_comut.add_categorical_data(mutect_data.loc[(mutect_data["Variant_Classification"].isin(set(step00.nonsynonymous_notations.values()))), ["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification"]].set_axis(labels=["sample", "category", "value"], axis="columns"), name="Mutation type", category_order=driver_data["Gene"], priority=["Frameshift indel"], mapping=step00.nonsynonymous_coloring)
    my_comut.add_bar_data(patient_data[["Tumor_Sample_Barcode", "TMB"]].set_axis(labels=["sample", "Counts"], axis="columns"), name="Mutation count", ylabel="TMB", mapping={"Counts": "purple"})
    # my_comut.add_side_bar_data(driver_data[["Gene", "-log10(P)"]].set_axis(labels=["category", "value"], axis="columns"), name="P-value", xlabel="-log10(P)", paired_name="Mutation type", position="left", mapping={"value": "olive"})
    my_comut.add_side_bar_data(driver_data[["Gene"] + stage_list].set_axis(labels=["category"] + stage_list, axis="columns"), name="Mutation rate", xlabel="Present Rate", paired_name="Mutation type", position="left", stacked=True, mapping=step00.stage_color_code)

    my_comut.plot_comut(x_padding=0.04, y_padding=0.04, tri_padding=0.03, figsize=(len(args.input) + 10, driver_data.shape[0] * 2))
    my_comut.add_unified_legend()
    my_comut.figure.savefig(args.output, bbox_inches="tight")
