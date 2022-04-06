"""
aggregate_mutect_sequenza.py: aggregate mutect MAF files with Sequenza CNV data
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

watching = "depth.ratio"
sequenza_data = pandas.DataFrame()


def get_data(file_name: str) -> pandas.DataFrame:
    data = pandas.read_csv(file_name, sep="\t", usecols=["chromosome", "start.pos", "end.pos", watching]).dropna(axis="index")
    data["sample"] = file_name.split("/")[-2]
    return data


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


def query(sample: str, chromosome: str, start: int, end: int) -> float:
    a = list()
    weights = list()

    tmp_data = sequenza_data.loc[(sequenza_data["sample"] == sample) & (sequenza_data["chromosome"] == chromosome) & (sequenza_data["start.pos"] <= start) & (end <= sequenza_data["end.pos"]), :]
    for index, row in tmp_data.iterrows():
        a.append(row["depth.ratio"])
        weights.append(end - start + 1)

    tmp_data = sequenza_data.loc[(sequenza_data["sample"] == sample) & (sequenza_data["chromosome"] == chromosome) & (start <= sequenza_data["end.pos"]) & (sequenza_data["end.pos"] <= end), :]
    for index, row in tmp_data.iterrows():
        a.append(row["depth.ratio"])
        weights.append(row["end.pos"] - start + 1)

    tmp_data = sequenza_data.loc[(sequenza_data["sample"] == sample) & (sequenza_data["chromosome"] == chromosome) & (start <= sequenza_data["start.pos"]) & (sequenza_data["start.pos"] <= end), :]
    for index, row in tmp_data.iterrows():
        a.append(row["depth.ratio"])
        weights.append(end - row["start.pos"] + 1)

    tmp_data = sequenza_data.loc[(sequenza_data["sample"] == sample) & (sequenza_data["chromosome"] == chromosome) & (start <= sequenza_data["start.pos"]) & (sequenza_data["end.pos"] <= end), :]
    for index, row in tmp_data.iterrows():
        a.append(row["depth.ratio"])
        weights.append(row["end.pos"] - row["start.pos"] + 1)

    if (not a) and (not weights):
        return 1.0
    else:
        return numpy.average(a, weights=weights)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("driver", help="MutEnricher Fisher enrichment output", type=str)
    parser.add_argument("census", help="Cancer gene census CSV file", type=str)
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output file", type=str)
    parser.add_argument("--cnv", help="CNV result(s)", type=str, nargs="+", required=True)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.05)
    parser.add_argument("--ratio", help="Ratio threshold for CNV", type=float, default=0.2)

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
    elif list(filter(lambda x: not x.endswith("sample_segments.txt"), args.cnv)):
        raise ValueError("CNV must end with sample_segments.txt!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.p < 1):
        raise ValueError("P-values must be (0, 1)")
    elif not (0 < args.ratio < 1):
        raise ValueError("ratio must be (0, 1)")

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
        mutect_data["Tumor_Sample_Barcode"] = pool.map(step00.get_id, mutect_data["Tumor_Sample_Barcode"])
    mutect_data["Variant_Classification"] = list(map(lambda x: step00.nonsynonymous_notations[x] if (x in step00.nonsynonymous_notations) else "Synonymous", mutect_data["Variant_Classification"]))
    mutect_data["Cancer_Stage"] = list(map(step00.get_long_sample_type, mutect_data["Tumor_Sample_Barcode"]))
    print(mutect_data)

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
    print(driver_data)

    with multiprocessing.Pool(args.cpus) as pool:
        sequenza_data = pandas.concat(objs=pool.map(get_data, args.cnv), axis="index", copy=False, ignore_index=True, verify_integrity=True)
    print(sequenza_data)

    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data["CNV"] = pool.starmap(query, mutect_data[["Tumor_Sample_Barcode", "Chromosome", "Start_Position", "End_Position"]].itertuples(index=False, name=None))
    print(min(mutect_data["CNV"]), numpy.mean(mutect_data["CNV"]), max(mutect_data["CNV"]))

    loss_data = mutect_data.loc[(mutect_data["CNV"] <= (1 - args.ratio))]
    loss_data["Variant_Classification"] = "CNV loss"

    gain_data = mutect_data.loc[(mutect_data["CNV"] >= (1 + args.ratio))]
    gain_data["Variant_Classification"] = "CNV gain"

    cnv_data = pandas.concat([loss_data, gain_data], ignore_index=True)
    print(cnv_data)

    mutect_data = pandas.concat([mutect_data, cnv_data], ignore_index=True)
    print(mutect_data)

    stage_list = list(filter(lambda x: x in list(mutect_data["Cancer_Stage"]), reversed(step00.long_sample_type_list)))
    print(stage_list)

    for stage in tqdm.tqdm(stage_list):
        snv_counter: collections.Counter = collections.Counter(mutect_data.loc[mutect_data["Cancer_Stage"] == stage].drop_duplicates(subset=["Hugo_Symbol", "Tumor_Sample_Barcode"])["Hugo_Symbol"])
        driver_data[stage] = list(map(lambda x: snv_counter[x] / len(args.input), driver_data["Gene"]))
    driver_data["-log10(P)"] = -1 * numpy.log10(driver_data["Fisher_pval"])
    print(driver_data)

    cnv_present_data = pandas.DataFrame()
    cnv_present_data["Hugo_Symbol"] = sorted(driver_data["Gene"])
    for stage in tqdm.tqdm(stage_list):
        cnv_counter: collections.Counter = collections.Counter(cnv_data.loc[cnv_data["Cancer_Stage"] == stage].drop_duplicates(subset=["Hugo_Symbol", "Tumor_Sample_Barcode"])["Hugo_Symbol"])
        cnv_present_data[stage] = list(map(lambda x: cnv_counter[x] / len(args.input), cnv_present_data["Hugo_Symbol"]))
    print(cnv_present_data)

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
    patient_data["TMB"] = list(map(lambda x: mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == x)].shape[0] / step00.WES_length * 10 ** 6, my_comut.samples))
    print(patient_data)

    if args.patient:
        my_comut.add_sample_indicators(patient_data[["Tumor_Sample_Barcode", "Patient"]].set_axis(labels=["sample", "group"], axis="columns"), name="Same patient")

    my_comut.add_categorical_data(patient_data[["Tumor_Sample_Barcode", "Collection_Type_category", "Collection_Type_value"]].set_axis(labels=["sample", "category", "value"], axis="columns"), name="Stage", mapping=step00.stage_color_code)
    my_comut.add_categorical_data(mutect_data.loc[(mutect_data["Variant_Classification"].isin(set(step00.nonsynonymous_notations.values()) | {"CNV loss", "CNV gain"})), ["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification"]].set_axis(labels=["sample", "category", "value"], axis="columns"), name="Mutation type", category_order=driver_data["Gene"], priority=["Frameshift indel"], mapping=step00.cnv_coloring, borders={"CNV loss", "CNV gain"})
    my_comut.add_bar_data(patient_data[["Tumor_Sample_Barcode", "TMB"]].set_axis(labels=["sample", "Counts"], axis="columns"), name="Mutation count", ylabel="TMB", mapping={"Counts": "purple"})
    # my_comut.add_side_bar_data(driver_data[["Gene"] + stage_list].set_axis(labels=["category"] + stage_list, axis="columns"), name="Mutation rate", xlabel="SNV Present Rate", paired_name="Mutation type", position="left", stacked=True, mapping=step00.stage_color_code)
    my_comut.add_side_bar_data(cnv_present_data.set_axis(labels=["category"] + stage_list, axis="columns"), name="Mutation rate", xlabel="CNV Present Rate", paired_name="Mutation type", position="left", stacked=True, mapping=step00.stage_color_code)

    my_comut.plot_comut(x_padding=0.04, y_padding=0.04, tri_padding=0.03, figsize=(len(args.input) + 5, driver_data.shape[0] * 2))
    my_comut.add_unified_legend()
    my_comut.figure.savefig(args.output, bbox_inches="tight")
