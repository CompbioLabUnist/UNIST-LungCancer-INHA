"""
calculate_SharedMutation_gene.py: calculate the number of of shared mutations in bar plot
"""
import argparse
import collections
import multiprocessing
import tarfile
import typing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Shared mutation information TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data w/ Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--percentage", help="Percentage of patients to include", type=float, default=0.25)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.05)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0.0 < args.percentage < 0.5):
        raise ValueError("Percentage must be (0.0, 0.5)!!")
    elif not (0 < args.p < 1):
        raise ValueError("P-values must be (0, 1)")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(len(patients), sorted(patients))

    input_data = input_data.loc[(input_data["Patient"].isin(patients))]
    print(input_data)

    filtered_data = input_data.loc[(input_data[step00.nonsynonymous_column].isin(step00.nonsynonymous_mutations))]
    print(filtered_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    for MSP in tqdm.tqdm(step00.sharing_columns[:1]):
        clinical_data = clinical_data.sort_values(MSP)

        lower_bound, higher_bound = numpy.quantile(clinical_data[MSP], args.percentage), numpy.quantile(clinical_data[MSP], 1.0 - args.percentage)

        lower_precancer_list = list(clinical_data.loc[(clinical_data[MSP] < lower_bound), f"{MSP}-sample"])
        higher_precancer_list = list(clinical_data.loc[(clinical_data[MSP] > higher_bound), f"{MSP}-sample"])

        lower_patients = list(map(step00.get_patient, lower_precancer_list))
        higher_patients = list(map(step00.get_patient, higher_precancer_list))

        lower_patients_length = len(lower_patients)
        higher_patients_length = len(higher_patients)

        if "SYN" in MSP:
            other_counter = collections.Counter(input_data.loc[~(input_data["Patient"].isin(lower_patients + higher_patients)), "Hugo_Symbol"])
            data = input_data.loc[(input_data["Patient"].isin(lower_patients + higher_patients))]
        else:
            other_counter = collections.Counter(filtered_data.loc[~(filtered_data["Patient"].isin(lower_patients + higher_patients)), "Hugo_Symbol"])
            data = filtered_data.loc[(filtered_data["Patient"].isin(lower_patients + higher_patients))]

        lower_counter = collections.Counter(data.loc[(data["Patient"].isin(lower_patients)), "Hugo_Symbol"])
        higher_counter = collections.Counter(data.loc[(data["Patient"].isin(higher_patients)), "Hugo_Symbol"])

        print(lower_counter.most_common(10))
        print(higher_counter.most_common(10))
        print(other_counter.most_common(10))

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(figures):
            tar.add(f, arcname=f)
