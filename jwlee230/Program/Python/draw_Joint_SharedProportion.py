"""
draw_Joint_SharedProportion.py: draw Joint plot with Shared mutation proportion
"""
import argparse
import tarfile
import matplotlib
import matplotlib.pyplot
import scipy.stats
import seaborn
import pandas
import tqdm
import step00

survival_columns = ["Recurrence-Free Survival", "Overall Survival"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--column", help="Column for Mutation Shared Proportion", choices=step00.sharing_columns, default=step00.sharing_columns[0])

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .tsv!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    clinical_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(len(patients))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()

    for column in tqdm.tqdm(survival_columns):
        clinical_data[column] = list(map(int, clinical_data[column]))

        r, p = scipy.stats.pearsonr(clinical_data[args.column], clinical_data[column])

        g = seaborn.jointplot(data=clinical_data, x=args.column, y=column, kind="reg", height=24, dropna=True)
        g.fig.text(0.5, 0.75, "r={0:.3f}, p={1:.3f}".format(r, p), color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")

        figures.append(f"{column}.pdf")
        g.savefig(figures[-1])

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
