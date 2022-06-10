"""
aggregate_PathSeq_clinical.py: Aggregate PathSeq results with clinical separation
"""
import argparse
import itertools
import multiprocessing
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import numpy
import pandas
import tqdm
import step00

input_data = pandas.DataFrame()


def get_data(filename: str) -> pandas.DataFrame:
    data = pandas.read_csv(filename, sep="\t")
    data["ID"] = step00.get_id(filename)
    return data


def get_real_taxonomy(taxon: str) -> str:
    return taxon.split("|")[-1].replace("_", " ")


def query(sample: str, taxon: str) -> float:
    data = input_data.loc[(input_data["real_taxonomy"] == taxon) & (input_data["ID"] == sample), "score_normalized"]
    if data.empty:
        return 0.0
    else:
        return data.to_numpy()[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PathSeq results TSV file(s)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--level", choices=step00.PathSeq_type_list, type=str, required=True)
    parser.add_argument("--compare", help="Comparison grouping (type, control, case, ...)", type=str, nargs="+", default=["Smoking-Detail", "Never", "Ex", "Current"])
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_sorting = parser.add_mutually_exclusive_group(required=True)
    group_sorting.add_argument("--patient", help="Sorting by patient first", action="store_true", default=False)
    group_sorting.add_argument("--type", help="Sorting by type first", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("INPUT must end with .TSV!!")
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
    print(patients)

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))

    if args.patient:
        args.input.sort(key=step00.sorting)
    elif args.type:
        args.input.sort(key=step00.sorting_by_type)
    else:
        raise Exception("Something went wrong!!")

    sample_list = list(map(step00.get_id, args.input))
    print(len(sample_list), sample_list)

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="index", copy=False, ignore_index=True, verify_integrity=True)
        input_data["real_taxonomy"] = pool.map(get_real_taxonomy, input_data["taxonomy"])
    input_data = input_data.loc[(input_data["kingdom"] == "Bacteria") & (input_data["type"] == args.level)]
    print(input_data)

    taxa_list = sorted(set(input_data["real_taxonomy"]))
    taxa_coloring = dict(zip(taxa_list, itertools.cycle(matplotlib.colors.XKCD_COLORS)))

    output_data = pandas.DataFrame(index=sample_list, columns=taxa_list, dtype=float)
    with multiprocessing.Pool(args.cpus) as pool:
        for index in tqdm.tqdm(sample_list):
            output_data.loc[index, :] = pool.starmap(query, [(index, taxon) for taxon in taxa_list])

    for index in tqdm.tqdm(sample_list):
        output_data.loc[index, :] = output_data.loc[index, :] / sum(output_data.loc[index, :]) * 100

    taxa_list = sorted(taxa_list, key=lambda x: numpy.mean(output_data[x]), reverse=True)
    output_data = output_data[taxa_list]
    output_data["Subtype"] = list(map(step00.get_long_sample_type, list(output_data.index)))
    output_data[args.compare[0]] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), args.compare[0]], list(output_data.index)))
    print(output_data)

    order = list(filter(lambda x: x in set(output_data["Subtype"]), step00.long_sample_type_list))
    print(order)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, axs = matplotlib.pyplot.subplots(nrows=len(args.compare) - 1, ncols=len(order), sharey="row", figsize=(len(sample_list) / 3, 9 * (len(args.compare) - 1)), gridspec_kw={"width_ratios": list(map(lambda x: len(output_data.loc[(output_data["Subtype"] == x)]), order))})

    for i, compare in enumerate(args.compare[1:]):
        for j, subtype in enumerate(order):
            drawing_data = output_data.loc[(output_data[args.compare[0]] == compare) & (output_data["Subtype"] == subtype)]
            taxa_list.sort(key=lambda x: numpy.mean(drawing_data.loc[:, x]), reverse=True)
            drawing_data = drawing_data.loc[:, taxa_list]
            print(compare, subtype, taxa_list[:10])
            for k, taxon in enumerate(tqdm.tqdm(taxa_list)):
                labeling = (i == 0) and (j == 0) and (k < 5)
                if labeling:
                    axs[i][j].bar(range(drawing_data.shape[0]), height=drawing_data.iloc[:, j], bottom=numpy.sum(drawing_data.iloc[:, :j], axis=1), color=taxa_coloring[taxon], label=taxon)
                else:
                    axs[i][j].bar(range(drawing_data.shape[0]), height=drawing_data.iloc[:, j], bottom=numpy.sum(drawing_data.iloc[:, :j], axis=1), color=taxa_coloring[taxon])

            axs[i][j].set_ylim([0, 100])
            axs[i][j].set_xticks([])
            axs[i][j].grid(True)
            axs[i][j].set_xlabel(f"{drawing_data.shape[0]} {subtype}", fontsize="xx-small")
            axs[i][j].set_title(f"{args.compare[0]}: {compare}")

            if i == 0:
                axs[i].legend(loc="lower left", fontsize="xx-small")
                axs[i].set_ylabel(f"Proportion in {args.level} (%)")

    matplotlib.pyplot.tight_layout()
    fig.savefig(args.output)
    fig.savefig(args.output.replace(".pdf", ".png"))
    matplotlib.pyplot.close(fig)
