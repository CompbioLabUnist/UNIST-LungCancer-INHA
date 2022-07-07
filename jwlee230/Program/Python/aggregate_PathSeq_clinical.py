"""
aggregate_PathSeq_clinical.py: Aggregate PathSeq results with clinical separation
"""
import argparse
import itertools
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import numpy
import pandas
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PathSeq results TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--level", choices=step00.PathSeq_type_list, type=str, required=True)
    parser.add_argument("--compare", help="Comparison grouping (type, control, case, ...)", type=str, nargs="+", default=["Smoking-Detail", "Never", "Ex", "Current"])

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_sorting = parser.add_mutually_exclusive_group(required=True)
    group_sorting.add_argument("--patient", help="Sorting by patient first", action="store_true", default=False)
    group_sorting.add_argument("--type", help="Sorting by type first", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    output_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    output_data[args.compare[0]] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), args.compare[0]], list(output_data.index)))
    index = list(filter(lambda x: step00.get_patient(x) in patients, list(output_data.index)))

    if args.patient:
        index.sort(key=step00.sorting)
    elif args.type:
        index.sort(key=step00.sorting_by_type)
    else:
        raise Exception("Something went wrong!!")

    output_data = output_data.loc[index, :]
    print(output_data)

    taxa_list = list(output_data.columns)[:-2]
    taxa_coloring = dict(zip(taxa_list, itertools.cycle(matplotlib.colors.XKCD_COLORS)))

    order = list(filter(lambda x: x in set(output_data["Subtype"]), step00.long_sample_type_list))
    print(order)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, axs = matplotlib.pyplot.subplots(ncols=len(args.compare) - 1, nrows=len(order), figsize=(16 * (len(args.compare) - 1), 9 * len(order)), sharey=True)

    for j, compare in enumerate(args.compare[1:]):
        for i, subtype in enumerate(order):
            drawing_data = output_data.loc[(output_data[args.compare[0]] == compare) & (output_data["Subtype"] == subtype)]
            taxa_list.sort(key=lambda x: numpy.mean(drawing_data.loc[:, x]), reverse=True)
            samples = sorted(list(drawing_data.index), key=lambda x: tuple(drawing_data.loc[x, taxa_list].to_numpy()), reverse=True)
            drawing_data = drawing_data.loc[samples, taxa_list]
            print(compare, subtype, taxa_list[:10])
            for k, taxon in enumerate(tqdm.tqdm(taxa_list)):
                labeling = (k < 5)
                if labeling:
                    axs[i][j].bar(range(drawing_data.shape[0]), height=drawing_data.iloc[:, k], bottom=numpy.sum(drawing_data.iloc[:, :k], axis=1), color=taxa_coloring[taxon], label=taxon)
                else:
                    axs[i][j].bar(range(drawing_data.shape[0]), height=drawing_data.iloc[:, k], bottom=numpy.sum(drawing_data.iloc[:, :k], axis=1), color=taxa_coloring[taxon])

            axs[i][j].set_ylim([0, 100])
            axs[i][j].set_xticks([])
            axs[i][j].grid(True)
            axs[i][j].set_xlabel(f"{drawing_data.shape[0]} {subtype} Samples")
            axs[i][j].set_title(f"{args.compare[0]}: {compare}")
            if drawing_data.shape[0]:
                axs[i][j].legend(loc="lower left", fontsize="xx-small")

            if j == 0:
                axs[i][j].set_ylabel(f"Proportion in {args.level} (%)")

    matplotlib.pyplot.tight_layout()
    fig.savefig(args.output)
    fig.savefig(args.output.replace(".pdf", ".png"))
    matplotlib.pyplot.close(fig)
