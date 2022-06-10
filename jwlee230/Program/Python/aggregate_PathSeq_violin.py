"""
aggregate_PathSeq_violin.py: Aggregate PathSeq results as violin plot
"""
import argparse
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import statannotations.Annotator
import tqdm
import step00

input_data = pandas.DataFrame()
output_data = pandas.DataFrame()


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


def draw_violin(taxon: str) -> str:
    stage_order = list(filter(lambda x: (x in set(output_data["Stage"])) and (len(output_data.loc[(output_data["Stage"] == x)]) > 3), step00.long_sample_type_list))

    try:
        stat, p = scipy.stats.kruskal(*[output_data.loc[(output_data["Stage"] == stage), taxon] for stage in stage_order])
    except ValueError:
        return ""

    if p >= 0.05:
        return ""

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.violinplot(data=output_data, x="Stage", y=taxon, order=stage_order, palette=step00.stage_color_code, cut=1, linewidth=5, ax=ax)
    statannotations.Annotator.Annotator(ax, list(zip(stage_order, stage_order[1:])), data=output_data, x="Stage", y=taxon, order=stage_order, palette=step00.stage_color_code).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.ylabel(f"{taxon} (%)")
    matplotlib.pyplot.title(f"Kruskal-Wallis p={p:.3f}")
    matplotlib.pyplot.tight_layout()

    fig_name = taxon.replace(" ", "_") + ".pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PathSeq results TSV file(s)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
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
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
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
    sample_list = list(map(step00.get_id, args.input))
    print(len(sample_list), sample_list)

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="index", copy=False, ignore_index=True, verify_integrity=True)
        input_data["real_taxonomy"] = pool.map(get_real_taxonomy, input_data["taxonomy"])
    input_data = input_data.loc[(input_data["kingdom"] == "Bacteria") & (input_data["type"] == args.level)]
    print(input_data)

    taxa_list = sorted(set(input_data["real_taxonomy"]))

    output_data = pandas.DataFrame(index=sample_list, columns=taxa_list, dtype=float)
    with multiprocessing.Pool(args.cpus) as pool:
        for index in tqdm.tqdm(sample_list):
            output_data.loc[index, :] = pool.starmap(query, [(index, taxon) for taxon in taxa_list])

    for index in tqdm.tqdm(sample_list):
        output_data.loc[index, :] = output_data.loc[index, :] / sum(output_data.loc[index, :]) * 100

    output_data["Stage"] = list(map(step00.get_long_sample_type, list(output_data.index)))
    taxa_list = sorted(taxa_list, key=lambda x: numpy.mean(output_data[x]), reverse=True)
    print(output_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    with multiprocessing.Pool(args.cpus) as pool:
        figures = list(filter(None, pool.map(draw_violin, taxa_list[:1000])))

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
