"""
aggregate_PathSeq.py: Aggregate PathSeq results
"""
import argparse
import itertools
import multiprocessing
import typing
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

    input_dict: typing.Dict[str, typing.List[str]] = dict()

    if args.patient:
        input_dict = {x: [] for x in patients}
        for sample in tqdm.tqdm(sample_list):
            input_dict[step00.get_patient(sample)].append(sample)
    elif args.type:
        input_dict = {x: [] for x in step00.long_sample_type_list}
        for sample in tqdm.tqdm(sample_list):
            input_dict[step00.get_long_sample_type(sample)].append(sample)
    else:
        raise Exception("Something went wrong!!")

    for key in tqdm.tqdm(list(input_dict.keys())):
        if not input_dict[key]:
            del input_dict[key]

    print(list(map(len, input_dict)))
    print(input_dict)

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
    print(output_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, axs = matplotlib.pyplot.subplots(nrows=1, ncols=len(input_dict), sharey="row", figsize=(64, 18), gridspec_kw={"width_ratios": list(map(len, input_dict.values()))})

    for i, (key, samples) in enumerate(input_dict.items()):
        samples.sort(key=lambda x: tuple(output_data.loc[x, :].to_numpy()), reverse=True)
        taxa_list.sort(key=lambda x: numpy.mean(output_data.loc[samples, x]), reverse=True)
        for j, taxon in enumerate(tqdm.tqdm(taxa_list)):
            labeling = (i == 0) and (j < 5)
            if labeling:
                axs[i].bar(range(len(samples)), height=output_data.loc[samples, :].iloc[:, j], bottom=numpy.sum(output_data.loc[samples, :].iloc[:, :j], axis=1), color=taxa_coloring[taxon], label=taxon)
            else:
                axs[i].bar(range(len(samples)), height=output_data.loc[samples, :].iloc[:, j], bottom=numpy.sum(output_data.loc[samples, :].iloc[:, :j], axis=1), color=taxa_coloring[taxon])

        axs[i].set_ylim([0, 100])
        axs[i].set_xticks([])
        axs[i].grid(True)
        axs[i].set_xlabel(f"{len(samples)} {key}")

        if i == 0:
            axs[i].legend(loc="lower left", fontsize="xx-small")

        if i == 0:
            axs[i].set_ylabel(f"Proportion in {args.level} (%)")

    matplotlib.pyplot.tight_layout()
    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
