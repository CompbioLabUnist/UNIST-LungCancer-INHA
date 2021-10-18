"""
aggregate_sequenza_correlation.py: aggregate sequenza results by correlation
"""
import argparse
import multiprocessing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import step00

centromere_dict = dict()
clinical_columns = ["Age", "Diagnosis_Size_OP", "Recurrence-Free Survial", "Overall Survival"]


def cut_ratio(value: float) -> float:
    if 0 <= value:
        return value
    elif value < 0:
        return 0
    else:
        raise ValueError("Something went wrong!!")


def get_data(file_name: str) -> pandas.DataFrame:
    data = pandas.read_csv(file_name, sep="\t", usecols=["chromosome", "start.pos", "end.pos", "depth.ratio"]).dropna(axis="index")
    data["sample"] = file_name.split("/")[-2]
    data["depth.ratio"] = list(map(cut_ratio, data["depth.ratio"]))
    return data


def is_centromere(chromosome: str, start: int, end: int) -> str:
    if chromosome not in centromere_dict:
        raise ValueError("Invalid chromosome name!!")

    if end < centromere_dict[row["chrom"]]["start"]:
        return chromosome + ".p"
    elif start > centromere_dict[row["chrom"]]["end"]:
        return chromosome + ".q"
    elif centromere_dict[row["chrom"]]["start"] < start < end < centromere_dict[row["chrom"]]["end"]:
        return "centromere"
    elif start < centromere_dict[row["chrom"]]["start"] < centromere_dict[row["chrom"]]["end"] < end:
        return "centromere"
    elif start < centromere_dict[row["chrom"]]["start"] < end:
        return chromosome + ".p"
    elif start < centromere_dict[row["chrom"]]["end"] < end:
        return chromosome + ".q"
    else:
        return "centromere"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Sequenza output segments.txt file(s)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("size", help="SIZE file", type=str)
    parser.add_argument("centromere", help="Centromere file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)
    parser.add_argument("--threshold", help="Threshold for gain/loss", type=float, default=0.2)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".txt"), args.input)):
        raise ValueError("INPUT must end with .TXT!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.threshold < 1):
        raise ValueError("Threshold must be (0, 1)")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)
    print(sorted(clinical_data.columns))

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    args.input = list(filter(lambda x: step00.get_patient(x.split("/")[-2]) in patients, args.input))

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="index", copy=False, ignore_index=True)
    input_data = input_data.loc[~(input_data["chromosome"] == "chrX")]
    print(input_data)

    chromosome_list = list(filter(lambda x: x in set(input_data["chromosome"]), step00.chromosome_list))
    sample_list = sorted(set(input_data["sample"]), key=step00.sorting)
    primary_cancer_list = list(filter(lambda x: step00.get_long_sample_type(x) == "Primary", sample_list))
    precancer_list = list(filter(lambda x: step00.get_long_sample_type(x) != "Primary", sample_list))
    print(chromosome_list)
    print(len(sample_list), sample_list)

    size_data = pandas.read_csv(args.size, sep="\t", header=None, names=["chromosome", "length"]).set_index(keys="chromosome", verify_integrity=True)
    print(size_data)

    centromere_data = pandas.read_csv(args.centromere, sep="\t")
    print(centromere_data)

    for index, row in centromere_data.iterrows():
        if row["chrom"] not in centromere_dict:
            centromere_dict[row["chrom"]] = {"start": row["chromStart"], "end": row["chromEnd"]}
        else:
            if row["chromStart"] < centromere_dict[row["chrom"]]["start"]:
                centromere_dict[row["chrom"]]["start"] = row["chromStart"]
            if row["chromEnd"] > centromere_dict[row["chrom"]]["end"]:
                centromere_dict[row["chrom"]]["end"] = row["chromEnd"]
    print(centromere_dict)

    with multiprocessing.Pool(args.cpus) as pool:
        input_data["detail_chromosome"] = pool.starmap(is_centromere, input_data[["chromosome", "start.pos", "end.pos"]].itertuples(index=False, name=None))
    input_data = input_data.loc[~(input_data["detail_chromosome"] == "centromere")]
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)
