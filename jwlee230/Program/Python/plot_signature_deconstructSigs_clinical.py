"""
plot_signature_deconstructSigs_clinical.py: violin plot cancer signature from deconstructSigs by stage with clinical data
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
compare_order: typing.List[typing.Tuple[typing.Tuple[str, str], typing.Tuple[str, str]]] = list()
hue_order: typing.List[str] = list()


def draw_violin(signature: str, clinical: str) -> str:
    try:
        stat, p = scipy.stats.kruskal(*[input_data.loc[(input_data["Stage"] == stage) & (input_data[clinical] == clinical_value), signature] for stage, clinical_value in itertools.product(order, hue_order)])
    except ValueError:
        _, p = 0.0, 1.0

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.violinplot(data=input_data, x=clinical, order=order, y=signature, hue="Stage", hue_order=hue_order, palette=step00.stage_color_code, inner="box", cut=1, linewidth=5, ax=ax)
    statannotations.Annotator.Annotator(ax, compare_order, data=input_data, x=clinical, order=order, y=signature, hue="Stage", hue_order=hue_order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.ylabel("Count")
    matplotlib.pyplot.title(f"{signature}: Kruskal-Wallis p={p:.3f}")
    matplotlib.pyplot.tight_layout()

    fig_name = f"{signature}.pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Signature TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--compare", help="Comparison grouping (type, control, case)", type=str, nargs="+", default=["Recurrence", "NO", "YES"])
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical data must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data.columns = list(map(lambda x: x.replace("Signature.", "SBS"), list(input_data.columns)))
    signatures = list(filter(lambda x: len(set(input_data[x])) > 1, list(input_data.columns)))
    print(input_data)

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC") & (clinical_data[args.compare[0]].isin(args.compare[1:]))].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC") & (clinical_data[args.compare[0]].isin(args.compare[1:]))].index)
    else:
        raise Exception("Something went wrong!!")

    input_data = input_data.loc[sorted(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index)), key=step00.sorting_by_type), :]
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    input_data[args.compare[0]] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), args.compare[0]], list(input_data.index)))
    print(input_data)

    order = args.compare[1:]
    hue_order = list(filter(lambda x: all([(not input_data.loc[(input_data["Stage"] == x) & (input_data[args.compare[0]] == c)].empty) for c in order]), step00.long_sample_type_list))
    compare_order = [((c, s1), (c, s2)) for s1, s2 in itertools.combinations(hue_order, r=2) for c in order] + [((c1, s), (c2, s)) for c1, c2 in itertools.combinations(order, r=2) for s in hue_order]
    print(order)
    print(hue_order)
    print(compare_order)

    with multiprocessing.Pool(args.cpus) as pool:
        figures = list(filter(None, pool.starmap(draw_violin, [(signature, args.compare[0]) for signature in signatures])))

    with tarfile.open(name=args.output, mode="w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
