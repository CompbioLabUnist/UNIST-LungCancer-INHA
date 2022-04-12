"""
plot_signature_bar.py: Plot signature from deconstructSigs in bar graph
"""
import argparse
import itertools
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import pandas
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Signature TSV file", type=str)
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_type = parser.add_mutually_exclusive_group(required=True)
    group_type.add_argument("--absolute", help="Absolute bar graph", action="store_true", default=False)
    group_type.add_argument("--relative", help="Relative bar graph", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical data must end with .csv!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data.columns = list(map(lambda x: x.replace("Signature.", "SBS"), list(input_data.columns)))
    signatures = list(input_data.columns)
    print(input_data)

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(sorted(patients))

    input_data = input_data.loc[sorted(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index)), key=step00.sorting_by_type), :]
    signatures = list(filter(lambda x: len(set(input_data[x])) > 1, list(input_data.columns)))
    input_data["Total"] = input_data.sum(axis="columns")
    print(input_data)
    print(signatures)

    if args.absolute:
        input_data.sort_values(by=sorted(signatures, key=lambda x: sum(input_data.loc[:, x]), reverse=True), ascending=False, inplace=True)
        input_data.sort_values(by="Total", kind="stable", ascending=False, inplace=True)
        input_data = input_data.loc[:, signatures]
    elif args.relative:
        for index in list(input_data.index):
            input_data.loc[index, :] = input_data.loc[index, :] / input_data.loc[index, "Total"]
        input_data.sort_values(by=sorted(signatures, key=lambda x: sum(input_data.loc[:, x]), reverse=True), ascending=False, inplace=True)
        input_data = input_data.loc[:, signatures]
    else:
        raise Exception("Something went wrong!!")
    input_data["Subtype"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

    for j, (column, color) in tqdm.tqdm(list(enumerate(zip(signatures, itertools.cycle(matplotlib.colors.TABLEAU_COLORS))))):
        matplotlib.pyplot.bar(range(input_data.shape[0]), list(input_data.loc[:, column]), bottom=input_data.iloc[:, :j].sum(axis="columns"), color=color, edgecolor=color, label=column)

    matplotlib.pyplot.xlabel("{0} Samples".format(input_data.shape[0]))
    if args.absolute:
        matplotlib.pyplot.ylabel("Counts")
    elif args.relative:
        matplotlib.pyplot.ylabel("Proportion")
    else:
        raise Exception("Something went wrong!!")
    matplotlib.pyplot.xticks([])
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.legend()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
