"""
draw_scatter_RSEM_MSP-ratio.py: draw scatter plots with DEG-ratio & MSP
"""
import argparse
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import tqdm
import step00

input_data = pandas.DataFrame()
expression_data = pandas.DataFrame()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input DEG-MSP TSV file", type=str)
    parser.add_argument("expression", help="Expression TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--r", help="r-value threshold", type=float, default=0.3)
    parser.add_argument("--slope", help="Slope threshold", type=float, default=5)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.expression.endswith(".tsv"):
        raise ValueError("Expression must end with .TSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif not (0 < args.r < 1):
        raise ValueError("r-value must be in (0, 1)!!")
    elif args.slope <= 0:
        raise ValueError("Slope must be positive!!")
    elif args.cpus < 1:
        raise ValueError("Number of CPUs must be positive!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)

    clinical_data: pandas.DataFrame = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    patients = set(clinical_data.index)
    print(clinical_data)

    expression_data = pandas.read_csv(args.expression, sep="\t", index_col=0).T
    expression_data = expression_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(expression_data.index))), :]
    expression_data["Stage"] = list(map(step00.get_long_sample_type, list(expression_data.index)))
    for column in tqdm.tqdm(step00.sharing_columns):
        expression_data[column] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), column], list(expression_data.index)))
    print(expression_data)

    stages = ["Precancer", "Primary"]

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    with multiprocessing.Pool(processes=args.cpus) as pool:
        for MSP in tqdm.tqdm(step00.sharing_columns[:1]):
            if f"{MSP}-log10(abs(slope))" not in set(input_data.columns):
                continue
