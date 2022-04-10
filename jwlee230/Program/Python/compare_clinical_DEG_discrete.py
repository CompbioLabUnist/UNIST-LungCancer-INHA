"""
compare_clinical_DEG_discrete.py: compare discrete clinical value with DEG
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

input_data = pandas.DataFrame()
column = ""
column_order: typing.List[str] = list()
hue_order: typing.List[str] = list()
p_threshold = 0.05


def run_1(gene: str) -> str:
    for c1, c2 in itertools.combinations(column_order, r=2):
        stat, p = scipy.stats.mannwhitneyu(input_data.loc[(input_data[column] == c1), gene], input_data.loc[(input_data[column] == c2), gene])
        if p < p_threshold:
            break
    else:
        return ""

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.violinplot(data=input_data, x=column, y=gene, order=column_order)
    statannotations.Annotator.Annotator(ax, list(itertools.combinations(column_order, 2)), data=input_data, x=column, y=gene, order=column_order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.xlabel(f"{column}")
    matplotlib.pyplot.ylabel(f"{gene}")
    matplotlib.pyplot.tight_layout()

    fig_name = f"{gene}_All.pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


def run_2(gene: str) -> str:
    flag = True
    for stage in hue_order:
        if not flag:
            break

        for c1, c2 in list(itertools.combinations(column_order, 2)):
            stat, p = scipy.stats.mannwhitneyu(input_data.loc[(input_data[column] == c1) & (input_data["Stage"] == stage), gene], input_data.loc[(input_data[column] == c2) & (input_data["Stage"] == stage), gene])
            if p < p_threshold:
                flag = False
                break

    if flag:
        return ""

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.violinplot(data=input_data, x=column, y=gene, order=column_order, hue="Stage", hue_order=hue_order, palette=step00.stage_color_code)
    statannotations.Annotator.Annotator(ax, [((c1, stage), (c2, stage)) for c1, c2 in list(itertools.combinations(column_order, 2)) for stage in hue_order], data=input_data, x=column, y=gene, order=column_order, hue="Stage", hue_order=hue_order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.xlabel(f"{column}")
    matplotlib.pyplot.ylabel(f"{gene}")
    matplotlib.pyplot.tight_layout()

    fig_name = f"{gene}_Stage.pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input gene expression TSV files", type=str)
    parser.add_argument("clinical", help="Clinical data data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--column", help="Which clinical value to use?", type=str, required=True)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)
    parser.add_argument("--p", help="P-value threshold", type=float, default=1e-2)

    group_histology = parser.add_mutually_exclusive_group(required=True)
    group_histology.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_histology.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .TSV!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output file must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("Number of CPUs must be positive!!")
    elif not (0 < args.p < 1):
        raise ValueError("P-value must be (0, 1)!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0).T
    patients = sorted(input_data.index, key=step00.sorting_by_type)
    genes = list(input_data.columns)
    print(input_data)

    clinical_data = step00.get_clinical_data(args.clinical)
    print(sorted(clinical_data.columns))
    assert args.column in set(clinical_data.columns)
    print(clinical_data)

    if args.SQC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
        patients = list(filter(lambda x: step00.get_patient(x) in histology, patients))
    elif args.ADC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
        patients = list(filter(lambda x: step00.get_patient(x) in histology, patients))
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    input_data = input_data.loc[patients, :]
    input_data[args.column] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), args.column], patients))
    input_data["Stage"] = list(map(step00.get_long_sample_type, patients))
    print(input_data)

    p_threshold = args.p
    column = args.column
    column_order = sorted(set(clinical_data[column]))
    hue_order = step00.long_sample_type_list[:]
    for order in column_order:
        hue_order = list(filter(lambda x: x in list(input_data.loc[(input_data[column] == order), "Stage"]), hue_order))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    with multiprocessing.Pool(args.cpus) as pool:
        figures = pool.map(run_1, genes)
        figures += pool.map(run_2, genes)

    figures = list(filter(None, figures))

    with tarfile.open(args.output, 'w') as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
