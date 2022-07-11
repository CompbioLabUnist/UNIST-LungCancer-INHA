"""
aggregate_PathSeq_violin_MutationSharedProportion.py: Aggregate PathSeq results as violin plot with Mutation Shared Proportion
"""
import argparse
import itertools
import multiprocessing
import tarfile
import typing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import statannotations.Annotator
import tqdm
import step00

output_data = pandas.DataFrame()
compare = ["Mutation Shared Proportion (Lower/Higher)", "Lower", "Higher"]


def draw_violin(taxon: str) -> str:
    stage_order: typing.List[str] = list(filter(lambda x: all([not output_data.loc[(output_data["Subtype"] == x) & (output_data[compare[0]] == comparing)].empty for comparing in compare[1:]]), step00.long_sample_type_list))

    try:
        stat, p = scipy.stats.kruskal(*[output_data.loc[(output_data["Subtype"] == stage) & (output_data[compare[0]] == comparing), taxon] for stage, comparing in itertools.product(stage_order, compare[1:])])
    except ValueError:
        return ""

    if p >= 0.01:
        return ""

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.violinplot(data=output_data, x="Subtype", y=taxon, order=stage_order, hue=compare[0], hue_order=compare[1:], cut=1, linewidth=5, ax=ax)
    try:
        statannotations.Annotator.Annotator(ax, list(map(lambda x: ((x[0], x[1][0]), (x[0], x[1][1])), itertools.product(stage_order, itertools.combinations(compare[1:], r=2)))), data=output_data, x="Subtype", y=taxon, order=stage_order, hue=compare[0], hue_order=compare[1:]).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()
    except Exception:
        pass

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
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--column", help="Column for Mutation Shared Proportion", choices=step00.sharing_columns, default=step00.sharing_columns[0])
    parser.add_argument("--level", choices=step00.PathSeq_type_list, type=str, required=True)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

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
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

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
        threshold = numpy.median(clinical_data[args.column])
    elif args.mean:
        threshold = numpy.mean(clinical_data[args.column])
    else:
        raise Exception("Something went wrong!!")
    print(f"{threshold:.3f}")

    clinical_data[compare[0]] = list(map(lambda x: "Higher" if (x > threshold) else "Lower", clinical_data[args.column]))
    print(clinical_data)

    output_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    output_data = output_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(output_data.index))), :]
    output_data[compare[0]] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), compare[0]], list(output_data.index)))
    print(output_data)

    taxa_list = list(output_data.columns)[:-2]

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    with multiprocessing.Pool(args.cpus) as pool:
        figures = list(filter(None, pool.map(draw_violin, taxa_list)))
    print(len(figures))

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
