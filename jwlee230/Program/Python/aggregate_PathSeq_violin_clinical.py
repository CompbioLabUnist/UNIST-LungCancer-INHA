"""
aggregate_PathSeq_violin_clinical.py: Aggregate PathSeq results as violin plot with clinical data
"""
import argparse
import itertools
import multiprocessing
import tarfile
import typing
import matplotlib
import matplotlib.pyplot
import pandas
import scipy.stats
import seaborn
import statannotations.Annotator
import tqdm
import step00

output_data = pandas.DataFrame()
compare: typing.List[str] = list()


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
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--level", choices=step00.PathSeq_type_list, type=str, required=True)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)
    parser.add_argument("--compare", help="Comparison grouping (type, control, case, ...)", type=str, nargs="+", default=["Smoking-Detail", "Never", "Ex", "Current"])

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    comapre = args.compare

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
    output_data = output_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(output_data.index))), :]
    output_data[args.compare[0]] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), args.compare[0]], list(output_data.index)))
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
