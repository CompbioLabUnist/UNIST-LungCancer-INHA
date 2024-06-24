"""
draw_Joint_SharedProportion.py: draw Joint plot with Shared mutation proportion
"""
import argparse
import itertools
import tarfile
import matplotlib
import matplotlib.pyplot
import scipy.stats
import seaborn
import pandas
import tqdm
import step00

survivals = ["Recurrence-Free Survival", "Overall Survival", "Pack-Year"]


def get_middle(values):
    return (min(values) + max(values)) / 2


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

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
    print("Patients:", len(patients))
    print(sorted(clinical_data.columns))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    for MSP, survival in tqdm.tqdm(list(itertools.product(step00.sharing_columns, survivals))):
        clinical_data[survival] = list(map(int, clinical_data[survival]))

        r, p = scipy.stats.spearmanr(clinical_data[MSP], clinical_data[survival])

        g = seaborn.jointplot(data=clinical_data, x=MSP, y=survival, kind="reg", height=24, dropna=True)
        g.set_axis_labels(MSP, f"{survival} (days)")
        g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})

        figures.append(f"Joint-{MSP}-{survival}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.regplot(data=clinical_data, x=MSP, y=survival, fit_reg=True, scatter=True, ax=ax)

        matplotlib.pyplot.ylabel(survival)
        matplotlib.pyplot.text(get_middle(clinical_data[MSP]), get_middle(clinical_data[survival]), f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        matplotlib.pyplot.tight_layout()

        figures.append(f"Scatter-{MSP}-{survival}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
