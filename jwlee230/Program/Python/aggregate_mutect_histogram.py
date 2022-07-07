"""
aggregate_mutect_histogram.py: aggregate mutect MAF files as histogram
"""
import argparse
import itertools
import multiprocessing
import matplotlib
import matplotlib.pyplot
import seaborn
import pandas
import tqdm
import step00

query_data = dict()


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


def query(sample_a: str, sample_b: str) -> float:
    return len(query_data[sample_a] & query_data[sample_b]) / len(query_data[sample_a])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

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
    sample_list = list(map(step00.get_id, args.input))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
        mutect_data["Tumor_Sample_Barcode"] = pool.map(step00.get_id, mutect_data["Tumor_Sample_Barcode"])
    print(mutect_data)

    for sample in tqdm.tqdm(sample_list):
        query_data[sample] = set(mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == sample), ["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]].itertuples(index=False, name=None))

    with multiprocessing.Pool(args.cpus) as pool:
        rate_list = pool.starmap(query, itertools.product(sample_list, repeat=2))

    fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

    seaborn.histplot(data=rate_list, stat="probability", binwidth=0.01, kde=True, ax=ax)

    ax.figure.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
