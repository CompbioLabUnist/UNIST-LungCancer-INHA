"""
draw_distribution_SharedProportion.py: draw distribution of mutation shared proportion
"""
import argparse
import matplotlib
import matplotlib.pyplot
import numpy
import seaborn
import pandas
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutation Sharing Proportion input TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--column", help="Column for Mutation Shared Proportion", choices=step00.sharing_columns, default=step00.sharing_columns[0])

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .tsv!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    input_data: pandas.DataFrame = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)

    if args.SQC:
        input_data = input_data.loc[(input_data["Histology"] == "SQC")]
    elif args.ADC:
        input_data = input_data.loc[(input_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(input_data.index)
    print(patients)

    print("min:", min(input_data[args.column]))
    print("mean:", numpy.mean(input_data[args.column]))
    print("median:", numpy.median(input_data[args.column]))
    print("max:", max(input_data[args.column]))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    g = seaborn.displot(input_data, x=args.column, kind="hist", stat="probability", kde=True, height=18, aspect=16 / 9)

    g.tight_layout()
    g.savefig(args.output)
