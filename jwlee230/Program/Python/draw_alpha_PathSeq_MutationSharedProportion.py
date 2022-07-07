"""
draw_alpha_PathSeq_MutationSharedProportion.py: draw alpha-diversity by Mutation Shared Proportion
"""
import argparse
import itertools
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

compare = ["Mutation Shared Proportion", "Lower", "Higher"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PathSeq results TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

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

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data = input_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index))), :]
    input_data[compare[0]] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), compare[0]], list(input_data.index)))
    alphas = list(input_data.columns)[:-2]
    print(input_data)

    stage_order = list(filter(lambda x: all([not input_data.loc[(input_data["Subtype"] == x) & (input_data[compare[0]] == comparing)].empty for comparing in compare[1:]]), step00.long_sample_type_list))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()

    for alpha in tqdm.tqdm(alphas):
        try:
            stat, p = scipy.stats.kruskal(*[input_data.loc[(input_data["Subtype"] == stage) & (input_data[compare[0]] == comparing), alpha] for stage, comparing in itertools.product(stage_order, compare[1:])])
        except ValueError:
            p = 1.0

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        seaborn.violinplot(data=input_data, x="Subtype", y=alpha, order=stage_order, hue=compare[0], hue_order=compare[1:], cut=1, linewidth=5, ax=ax)
        try:
            statannotations.Annotator.Annotator(ax, list(map(lambda x: ((x[0], x[1][0]), (x[0], x[1][1])), itertools.product(stage_order, itertools.combinations(compare[1:], r=2)))), data=input_data, x="Subtype", y=alpha, order=stage_order, hue=compare[0], hue_order=compare[1:]).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()
        except Exception:
            pass

        matplotlib.pyplot.ylabel(f"{alpha.replace('_', ' ')}")
        matplotlib.pyplot.title(f"{compare[0]} (K.W. p={p:.3f})")
        matplotlib.pyplot.tight_layout()

        fig_name = f"{alpha}.pdf"
        fig.savefig(fig_name)
        matplotlib.pyplot.close(fig)

        figures.append(fig_name)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
