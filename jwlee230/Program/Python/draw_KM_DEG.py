"""
draw_KM_DEG.py: Draw survival plot with DEG
"""
import argparse
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import scipy.stats
import seaborn
import pandas
import tqdm
import step00

input_data = pandas.DataFrame()
clinical_data = pandas.DataFrame()
survival_columns = ["Recurrence-Free Survival", "Overall Survival"][1:]


def get_middle(values):
    return (min(values) + max(values)) / 2


def run_reg(gene: str, sharing: str, stage: str, survival: str) -> str:
    precancer_set = set(clinical_data[f"{sharing}-sample"])
    primary_set = set(map(step00.get_paired_primary, precancer_set))

    if stage == "Primary":
        drawing_data = input_data.loc[list(filter(lambda x: x in primary_set, list(input_data.index))), [gene]]
        color = "gray"
    elif stage == "Precancer":
        drawing_data = input_data.loc[list(filter(lambda x: x in precancer_set, list(input_data.index))), [gene]]
        color = "tab:pink"
    else:
        raise Exception("Something went wrong!!")

    if drawing_data.empty:
        return ""

    drawing_data[survival] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), survival], list(drawing_data.index)))

    r, p = scipy.stats.pearsonr(drawing_data[gene], drawing_data[survival])

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    seaborn.regplot(data=drawing_data, x=gene, y=survival, fit_reg=True, scatter=True, color=color, ax=ax)

    matplotlib.pyplot.text(get_middle(drawing_data[gene]), get_middle(drawing_data[survival]), f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
    matplotlib.pyplot.ylabel(f"{survival} (days)")
    matplotlib.pyplot.tight_layout()

    fig_name = f"reg-{gene}-{sharing}-{stage}-{survival}.pdf"
    matplotlib.pyplot.savefig(fig_name)
    matplotlib.pyplot.close()

    return fig_name


def run_lm(gene: str, sharing: str, stage: str, survival: str) -> str:
    precancer_set = set(clinical_data[f"{sharing}-sample"])
    primary_set = set(map(step00.get_paired_primary, precancer_set))

    if stage == "Primary":
        drawing_data = input_data.loc[list(filter(lambda x: x in primary_set, list(input_data.index))), [gene]]
        color = "gray"
    elif stage == "Precancer":
        drawing_data = input_data.loc[list(filter(lambda x: x in precancer_set, list(input_data.index))), [gene]]
        color = "tab:pink"
    else:
        raise Exception("Something went wrong!!")

    if drawing_data.empty:
        return ""

    drawing_data[survival] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), survival], list(drawing_data.index)))

    r, p = scipy.stats.pearsonr(drawing_data[gene], drawing_data[survival])

    g = seaborn.jointplot(data=drawing_data, x=gene, y=survival, kind="reg", height=24, ratio=5, dropna=True, color=color)
    g.set_axis_labels(gene, f"{survival} (days)")
    g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})

    fig_name = f"lm-{gene}-{sharing}-{stage}-{survival}.pdf"
    g.savefig(fig_name)
    matplotlib.pyplot.close(g.fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="DEG TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .tsv!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("Number of CPUs must be positive!!")

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(len(patients), patients)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0).T
    input_data = input_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index))), :]
    print(input_data)

    # gene_list = list(input_data.columns)
    gene_list = ["ELOB", "H2AC4", "KRI1", "MZT2B", "PRDX2", "PTGES3", "TJP2", "ZNF148", "MAPK8IP3"]
    print("Gene:", len(gene_list))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    with multiprocessing.Pool(processes=args.cpus) as pool:
        figures = list(filter(None, pool.starmap(run_reg, itertools.product(gene_list, step00.sharing_columns[:1], ["Precancer", "Primary"][:1], survival_columns))))

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
