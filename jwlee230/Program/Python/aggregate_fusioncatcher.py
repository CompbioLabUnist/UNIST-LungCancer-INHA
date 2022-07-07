"""
aggregate_fusioncatcher.py: Aggregate FusionCatcher results as heatmap plot
"""
import argparse
import multiprocessing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import step00


def read_maf(filename: str) -> pandas.DataFrame:
    data = pandas.read_csv(filename, sep="\t", low_memory=False)
    data["sample"] = step00.get_id(filename)
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="FusionCatcher output .tsv files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
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

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
    print(input_data)

    first_gene_set = set(input_data.iloc[:, 0])
    second_gene_set = set(input_data.iloc[:, 1])

    output_data = pandas.DataFrame(data=numpy.zeros((len(first_gene_set), len(second_gene_set))), index=sorted(first_gene_set), columns=sorted(second_gene_set), dtype=int)
    for index, row in input_data.iterrows():
        output_data.loc[row[0], row[1]] += 1
    print(output_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(len(first_gene_set), len(second_gene_set)))

    seaborn.heatmap(data=output_data, cmap="Reds", square=False, linecolor="k", linewidths=0.1, xticklabels=True, yticklabels=True, ax=ax)

    matplotlib.pyplot.xlabel("3' fusion gene")
    matplotlib.pyplot.ylabel("5' fusion gene")
    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
