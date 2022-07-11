"""
draw_KM_SharedProportion.py: Draw Kaplan-Meier plot with mutation shared proportion
"""
import argparse
import tarfile
import lifelines
import lifelines.statistics
import matplotlib
import matplotlib.pyplot
import numpy
import seaborn
import pandas
import tqdm
import step00

survival_columns = ["Recurrence-Free Survival", "Overall Survival"]
threshold = 365 * 5


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutation Sharing Proportion input TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--column", help="Column for Mutation Shared Proportion", choices=step00.sharing_columns, default=step00.sharing_columns[0])
    parser.add_argument("--cutting", help="Cutting follow-up up to 5 year", action="store_true", default=False)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_divide = parser.add_mutually_exclusive_group(required=True)
    group_divide.add_argument("--median", help="Divide as median", action="store_true", default=False)
    group_divide.add_argument("--mean", help="Divide as mean", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .tsv!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    input_data: pandas.DataFrame = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)

    if args.SQC:
        input_data = input_data.loc[(input_data["Histology"] == "SQC")]
    elif args.ADC:
        input_data = input_data.loc[(input_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(input_data.index)
    print(len(patients))

    if args.cutting:
        for column in tqdm.tqdm(survival_columns):
            input_data[column] = list(map(lambda x: threshold if (x > threshold) else x, input_data[column]))
    print(input_data)

    if args.median:
        threshold = numpy.median(input_data[args.column])
    elif args.mean:
        threshold = numpy.mean(input_data[args.column])
    else:
        raise Exception("Something went wrong!!")

    lower_data = input_data.loc[(input_data[args.column] <= threshold)]
    higher_data = input_data.loc[(input_data[args.column] > threshold)]
    print(lower_data)
    print(higher_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()

    for column in tqdm.tqdm(survival_columns):
        fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

        kmf = lifelines.KaplanMeierFitter()

        kmf.fit(lower_data[column], label=f"Lower Shared Proportion ({len(lower_data)} patients)")
        kmf.plot(ax=ax, ci_show=False, c="tab:blue")

        kmf.fit(higher_data[column], label=f"Higher Shared Proportion ({len(higher_data)} patients)")
        kmf.plot(ax=ax, ci_show=False, c="tab:red")

        p_value = lifelines.statistics.logrank_test(lower_data[column], higher_data[column]).p_value
        matplotlib.pyplot.text(max(input_data[column]) / 2, 1.0, f"p={p_value:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")

        matplotlib.pyplot.xlabel(f"{column} (Days)")
        matplotlib.pyplot.ylabel("Survival Rate")
        matplotlib.pyplot.title(f"Lower/Higher Threshold: {threshold:.3f}")
        matplotlib.pyplot.tight_layout()

        figures.append(f"{column.replace(' ', '-')}.pdf")
        matplotlib.pyplot.savefig(figures[-1])
        matplotlib.pyplot.close()

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
