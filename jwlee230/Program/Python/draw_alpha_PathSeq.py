"""
draw_alpha_PathSeq.py: draw alpha-diversity by stage
"""
import argparse
import tarfile
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

    parser.add_argument("input", help="PathSeq alpha-diversity results TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

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

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data = input_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index))), :]
    alphas = list(input_data.columns)[:-1]
    print(input_data)

    stage_order = list(filter(lambda x: (x in set(input_data["Subtype"])) and (len(input_data.loc[(input_data["Subtype"] == x)]) > 3), step00.long_sample_type_list))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()

    for alpha in tqdm.tqdm(alphas):
        try:
            stat, p = scipy.stats.kruskal(*[input_data.loc[(input_data["Subtype"] == stage), alpha] for stage in stage_order])
        except ValueError:
            p = 1.0

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        seaborn.violinplot(data=input_data, x="Subtype", y=alpha, order=stage_order, palette=step00.stage_color_code, cut=1, linewidth=5, ax=ax)
        try:
            statannotations.Annotator.Annotator(ax, list(zip(stage_order, stage_order[1:])), data=input_data, x="Subtype", y=alpha, order=stage_order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()
        except ValueError:
            pass

        matplotlib.pyplot.ylabel(f"{alpha.replace('_', ' ')}")
        matplotlib.pyplot.title(f"Kruskal-Wallis p={p:.3f}")
        matplotlib.pyplot.tight_layout()

        fig_name = f"{alpha}.pdf"
        fig.savefig(fig_name)
        matplotlib.pyplot.close(fig)

        figures.append(fig_name)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
