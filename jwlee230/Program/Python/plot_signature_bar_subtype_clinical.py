"""
plot_signature_bar_subtype.py: Plot signature in bar graph with subtype separation
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

    parser.add_argument("input", help="Signature TSV file (not necessarily TSV)", type=str)
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--compare", help="Comparison grouping (type, control, case)", type=str, nargs="+", default=["Smoking-Detail", "Never", "Ex", "Current"])

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_type = parser.add_mutually_exclusive_group(required=True)
    group_type.add_argument("--absolute", help="Absolute bar graph", action="store_true", default=False)
    group_type.add_argument("--relative", help="Relative bar graph", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".csv"):
        raise ValueError("Clinical data must end with .csv!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col="Samples")
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
    signatures = list(input_data.columns)
    input_data["Total"] = input_data.sum(axis="columns")
    print(input_data)
    print(signatures)

    if args.absolute:
        input_data.sort_values(by="Total", ascending=False, inplace=True)
        signatures.sort(key=lambda x: sum(input_data[x]), reverse=True)
    elif args.relative:
        for index in list(input_data.index):
            input_data.loc[index, :] = input_data.loc[index, :] / input_data.loc[index, "Total"]
        signatures.sort(key=lambda x: sum(input_data[x]), reverse=True)
        input_data.sort_values(by=signatures, ascending=False, inplace=True)
    else:
        raise Exception("Something went wrong!!")
    input_data = input_data.loc[:, signatures + ["Total"]]
    input_data["Subtype"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    input_data[args.compare[0]] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), args.compare[0]], list(input_data.index)))
    print(input_data)

    order = list(filter(lambda x: x in set(input_data["Subtype"]), step00.long_sample_type_list))
    print(order)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, axs = matplotlib.pyplot.subplots(ncols=len(args.compare) - 1, nrows=len(order), figsize=(16 * (len(args.compare) - 1), 9 * len(order)), sharey=True)

    for i, subtype in tqdm.tqdm(enumerate(order)):
        for j, clinical in enumerate(args.compare[1:]):
            drawing_data = input_data.loc[(input_data["Subtype"] == subtype) & (input_data[args.compare[0]] == clinical), :]

            if args.absolute:
                drawing_data.sort_values(by="Total", ascending=False, inplace=True)
            elif args.relative:
                drawing_data.sort_values(by=signatures, ascending=False, inplace=True)
            else:
                raise Exception("Something went wrong!!")

            for k, (column, color) in enumerate(zip(signatures, itertools.cycle(matplotlib.colors.TABLEAU_COLORS))):
                axs[i][j].bar(range(drawing_data.shape[0]), list(drawing_data.loc[:, column]), bottom=drawing_data.iloc[:, :k].sum(axis="columns"), color=color, edgecolor=color, label=column)

            axs[i][j].set_xlabel("{0} Samples".format(subtype))
            if args.absolute:
                axs[i][j].set_ylabel("Counts")
            elif args.relative:
                axs[i][j].set_ylabel("Proportion")
            else:
                raise Exception("Something went wrong!!")
            axs[i][j].set_xticks([])
            axs[i][j].grid(True)
            if i == 0 and j == len(args.compare) - 2:
                axs[i][j].legend()
            axs[i][j].set_title("{0}: {1}".format(args.compare[0], clinical))

    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
