"""
draw_Tcell_density_violin_deconvolution_clinical.py: draw T-cell density violin plot for deconvolution result with clinical values
"""
import argparse
import itertools
import matplotlib
import matplotlib.pyplot
import pandas
import scipy.stats
import seaborn
import statannotations.Annotator
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Deconvolution result TSV file (not necessarily TSV)", type=str)
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--column", help="Clinical data column name", nargs="+", default=["Recurrence", "NO", "YES"])

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        clincal_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clincal_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")

    patients = set(clinical_data.index)
    print(patients)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)

    input_data = input_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index))), step00.Tcell_list]
    input_data["T cell density"] = list(map(lambda x: sum(input_data.loc[x, step00.Tcell_list]), list(input_data.index)))
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    input_data[args.column[0]] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), args.column[0]], list(input_data.index)))
    print(input_data)

    for stage in tqdm.tqdm(set(input_data["Stage"])):
        if any([(len(input_data.loc[(input_data["Stage"] == stage) & (input_data[args.column[0]] == clinical)]) < 3) for clinical in args.column[1:]]):
            input_data = input_data.loc[~(input_data["Stage"] == stage)]
    print(input_data)

    stage_list = set(input_data["Stage"])
    order = list(filter(lambda x: x in stage_list, step00.long_sample_type_list))
    palette = list(map(lambda x: step00.stage_color_code[x], order))

    try:
        stat, p = scipy.stats.kruskal(*[input_data.loc[(input_data["Stage"] == stage) & (input_data[args.column[0]] == clinical), "T cell density"] for clinical in args.column[1:] for stage in order])
    except ValueError:
        stat, p = 0.0, 1.0

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.violinplot(data=input_data, x=args.column[0], y="T cell density", order=args.column[1:], hue="Stage", hue_order=order, palette=palette, cut=1, linewidth=5, ax=ax)
    statannotations.Annotator.Annotator(ax, [((clinical, a), (clinical, b)) for a, b in itertools.combinations(order, r=2) for clinical in args.column[1:]] + [((a, stage), (b, stage)) for a, b in zip(args.column[1:], args.column[2:]) for stage in order], data=input_data, x=args.column[0], y="T cell density", order=args.column[1:], hue="Stage", hue_order=order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.title(f"Kruskal-Wallis p={p:.3f}")
    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
