"""
aggregate_PathSeq_violin.py: Aggregate PathSeq results as violin plot
"""
import argparse
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import scipy.stats
import seaborn
import statannotations.Annotator
import tqdm
import step00

output_data = pandas.DataFrame()


def draw_violin(taxon: str) -> str:
    stage_order = list(filter(lambda x: (x in set(output_data["Subtype"])) and (len(output_data.loc[(output_data["Subtype"] == x)]) > 3), step00.long_sample_type_list))

    try:
        stat, p = scipy.stats.kruskal(*[output_data.loc[(output_data["Subtype"] == stage), taxon] for stage in stage_order])
    except ValueError:
        return ""

    if p >= 0.01:
        return ""

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.violinplot(data=output_data, x="Subtype", y=taxon, order=stage_order, palette=step00.stage_color_code, cut=1, linewidth=5, ax=ax)
    statannotations.Annotator.Annotator(ax, list(zip(stage_order, stage_order[1:])), data=output_data, x="Subtype", y=taxon, order=stage_order, palette=step00.stage_color_code).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.ylabel(f"{taxon} (%)")
    matplotlib.pyplot.title(f"Kruskal-Wallis p={p:.3f}")
    matplotlib.pyplot.tight_layout()

    fig_name = taxon.replace(" ", "_").replace("/", "_") + ".pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PathSeq results TSV file", type=str)
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

    if not args.input.endswith(".tsv"):
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

    output_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    index = list(filter(lambda x: step00.get_patient(x) in patients, list(output_data.index)))

    if args.patient:
        index.sort(key=step00.sorting)
    elif args.type:
        index.sort(key=step00.sorting_by_type)
    else:
        raise Exception("Something went wrong!!")

    output_data = output_data.loc[index, :]
    print(output_data)

    taxa_list = list(output_data.columns)[:-1]

    order = list(filter(lambda x: x in set(output_data["Subtype"]), step00.long_sample_type_list))
    print(order)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    with multiprocessing.Pool(args.cpus) as pool:
        figures = list(filter(None, pool.map(draw_violin, taxa_list)))

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
