"""
plot_signature_violin_relative.py: violin plot cancer signature by stage with relative count
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
order: typing.List[str] = list()


def draw_violin(signature: str) -> str:
    pairs = list()
    for stage1, stage2 in itertools.combinations(order, r=2):
        p = scipy.stats.mannwhitneyu(input_data.loc[(input_data["Stage"] == stage1), signature], input_data.loc[(input_data["Stage"] == stage2), signature])[1]
        if p < 0.05:
            pairs.append((stage1, stage2))

    try:
        stat, p = scipy.stats.kruskal(*[input_data.loc[(input_data["Stage"] == stage), signature] for stage in order])
    except ValueError:
        _, p = 0.0, 1.0

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    seaborn.violinplot(data=input_data, x="Stage", y=signature, order=order, palette=step00.stage_color_code, inner="box", linewidth=5, ax=ax)
    if pairs:
        statannotations.Annotator.Annotator(ax, pairs, data=input_data, x="Stage", y=signature, order=order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()

    matplotlib.pyplot.ylabel("Proportion")
    matplotlib.pyplot.title(f"{signature}: K.W. p={p:.3f}")
    matplotlib.pyplot.tight_layout()

    fig_name = f"{signature}.pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Signature TSV file (not necessarily TSV)", type=str)
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".csv"):
        raise ValueError("Clinical data must end with .csv!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = pandas.read_csv(args.input, sep="\t", index_col="Samples")
    signatures = list(input_data.columns)
    print(input_data)

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(sorted(patients))

    input_data = input_data.loc[sorted(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index)), key=step00.sorting_by_type), :]
    input_data["Total"] = input_data.sum(axis="columns")
    for index in tqdm.tqdm(list(input_data.index)):
        input_data.loc[index, :] = input_data.loc[index, :] / input_data.loc[index, "Total"]
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    order = list(filter(lambda x: x in set(input_data["Stage"]), step00.long_sample_type_list))
    print(input_data)
    print(order)

    with multiprocessing.Pool(args.cpus) as pool:
        figures = list(filter(None, pool.map(draw_violin, signatures)))

    with tarfile.open(name=args.output, mode="w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
