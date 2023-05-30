"""
aggregate_CNV_violin_ploidy.py: Violin plot of CNV data for cancer stages
"""
import argparse
import matplotlib
import matplotlib.pyplot
import pandas
import scipy.stats
import seaborn
import statannotations.Annotator
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="CNV sample.tsv file", type=str)
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

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
    input_data["Patient"] = list(map(step00.get_patient, list(input_data.index)))
    input_data = input_data.loc[(input_data["Patient"].isin(patients))]
    print(input_data)

    stage_list = list(filter(lambda x: x in set(input_data["Stage"]), step00.long_sample_type_list))
    print(stage_list)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    pairs = list()
    for stage_a, stage_b in zip(stage_list, stage_list[1:]):
        p = scipy.stats.mannwhitneyu(input_data.loc[(input_data["Stage"] == stage_a), "Ploidy"], input_data.loc[(input_data["Stage"] == stage_b), "Ploidy"])[1]
        if p < 0.05:
            pairs.append((stage_a, stage_b))

    try:
        stat, p = scipy.stats.kruskal(*[input_data.loc[(input_data["Stage"] == stage), "Ploidy"] for stage in stage_list])
    except ValueError:
        stat, p = 0.0, 1.0

    seaborn.violinplot(data=input_data, x="Stage", y="Ploidy", order=stage_list, inner="box", palette=step00.stage_color_code, cut=1, ax=ax)
    if pairs:
        statannotations.Annotator.Annotator(ax, pairs, data=input_data, x="Stage", y="Ploidy", order=stage_list).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()

    matplotlib.pyplot.title(f"Kruskal-Wallis p={p:.3f}")
    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
