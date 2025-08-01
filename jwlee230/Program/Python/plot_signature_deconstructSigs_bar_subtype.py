"""
plot_signature_deconstructSigs_bar_subtype.py: Plot signature from deconstructSigs in bar graph with subtype separation
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
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
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
    signatures = list(filter(lambda x: len(set(input_data[x])) > 1, list(input_data.columns)))
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
    input_data["Total"] = input_data.sum(axis="columns")
    print(input_data)
    print(signatures)

    input_data = input_data.loc[:, signatures + ["Total"]]
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    print(input_data)

    order = list(filter(lambda x: x in set(input_data["Stage"]), step00.long_sample_type_list))
    print(order)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, axs = matplotlib.pyplot.subplots(nrows=len(order), figsize=(32, 9 * len(order)), sharey=True)

    for i, subtype in tqdm.tqdm(list(enumerate(order))):
        drawing_data = input_data.loc[(input_data["Stage"] == subtype), signatures + ["Total"]].copy()

        if args.absolute:
            drawing_data.sort_values(by=sorted(signatures, key=lambda x: sum(drawing_data.loc[:, x]), reverse=True), ascending=False, inplace=True)
            drawing_data.sort_values(by="Total", kind="stable", ascending=False, inplace=True)
        elif args.relative:
            for index in list(drawing_data.index):
                drawing_data.loc[index, :] = drawing_data.loc[index, :] / drawing_data.loc[index, "Total"]
            drawing_data.sort_values(by=sorted(signatures, key=lambda x: sum(drawing_data.loc[:, x]), reverse=True), ascending=False, inplace=True)
        else:
            raise Exception("Something went wrong!!")

        for j, (column, color) in enumerate(zip(signatures, itertools.cycle(matplotlib.colors.XKCD_COLORS))):
            axs[i].bar(range(drawing_data.shape[0]), list(drawing_data.loc[:, column]), bottom=drawing_data.iloc[:, :j].sum(axis="columns"), color=color, edgecolor=color, label=column)

        if args.absolute:
            axs[i].set_ylabel("Counts")
        elif args.relative:
            axs[i].set_ylabel("Proportion")
        else:
            raise Exception("Something went wrong!!")

        if i == 0:
            axs[i].legend(ncol=3)

        axs[i].set_xlabel("{0} {1} Samples".format(drawing_data.shape[0], subtype))
        axs[i].set_xticks([])
        axs[i].grid(True)
        axs[i].set_title(subtype)

    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
