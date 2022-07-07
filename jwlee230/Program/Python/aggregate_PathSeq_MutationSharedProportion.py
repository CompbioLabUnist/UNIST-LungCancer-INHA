"""
aggregate_PathSeq_MutationSharedProportion.py: Aggregate PathSeq results with Mutation Shared Proportion
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

compare = ["Mutation Shared Proportion", "Lower", "Higher"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PathSeq results TSV file", type=str)
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--level", choices=step00.PathSeq_type_list, type=str, required=True)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_compare = parser.add_mutually_exclusive_group(required=True)
    group_compare.add_argument("--median", help="Compare median", action="store_true", default=False)
    group_compare.add_argument("--mean", help="Compare mean", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(patients)

    if args.median:
        threshold = numpy.median(clinical_data["Shared Proportion"])
    elif args.mean:
        threshold = numpy.mean(clinical_data["Shared Proportion"])
    else:
        raise Exception("Something went wrong!!")
    print(f"{threshold:.3f}")

    clinical_data[compare[0]] = list(map(lambda x: "Higher" if (x > threshold) else "Lower", clinical_data["Shared Proportion"]))
    print(clinical_data)

    output_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    output_data = output_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(output_data.index))), :]
    output_data[compare[0]] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), compare[0]], list(output_data.index)))
    print(output_data)

    taxa_list = list(output_data.columns)[:-2]
    taxa_coloring = dict(zip(taxa_list, itertools.cycle(matplotlib.colors.XKCD_COLORS)))

    order = list(filter(lambda x: x in set(output_data["Subtype"]), step00.long_sample_type_list))
    print(order)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, axs = matplotlib.pyplot.subplots(ncols=len(compare) - 1, nrows=len(order), figsize=(16 * (len(compare) - 1), 9 * len(order)), sharey=True)

    for j, comparing in enumerate(compare[1:]):
        for i, subtype in enumerate(order):
            drawing_data = output_data.loc[(output_data[compare[0]] == comparing) & (output_data["Subtype"] == subtype)]
            taxa_list.sort(key=lambda x: numpy.mean(drawing_data.loc[:, x]), reverse=True)
            samples = sorted(list(drawing_data.index), key=lambda x: tuple(drawing_data.loc[x, taxa_list].to_numpy()), reverse=True)
            drawing_data = drawing_data.loc[samples, taxa_list]
            print(comparing, subtype, taxa_list[:10])
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
            axs[i][j].set_title(f"{comparing}")
            if drawing_data.shape[0]:
                axs[i][j].legend(loc="lower left", fontsize="xx-small")

            if j == 0:
                axs[i][j].set_ylabel(f"Proportion in {args.level} (%)")

    matplotlib.pyplot.tight_layout()
    fig.savefig(args.output)
    fig.savefig(args.output.replace(".pdf", ".png"))
    matplotlib.pyplot.close(fig)
