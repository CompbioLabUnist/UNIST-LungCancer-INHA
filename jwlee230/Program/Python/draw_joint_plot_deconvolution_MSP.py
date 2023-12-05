"""
draw_joint_plot_deconvolution_MSP.py: draw joint plot from deconvolution results with MSP
"""
import argparse
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy
import seaborn
import tqdm
import step00

input_data = pandas.DataFrame()
replaced = " ()/"


def get_middle(values):
    return (max(values) + min(values)) / 2


def reg(stage, MSP, cell):
    if stage == "All":
        drawing_data = input_data.copy()
        color = "tab:blue"
    elif stage == "Precancer":
        drawing_data = input_data.loc[~(input_data["Stage"].isin({"Normal", "Primary"}))]
        color = "tab:pink"
    else:
        drawing_data = input_data.loc[(input_data["Stage"] == stage)]
        color = step00.stage_color_code[stage]

    if drawing_data.empty:
        return ""
    elif (numpy.std(drawing_data[MSP]) == 0) and (numpy.std(drawing_data[cell]) == 0):
        return ""

    slope, intercept, r, p, *_ = scipy.stats.linregress(drawing_data[MSP], drawing_data[cell])

    if p >= 0.05:
        return ""

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    seaborn.regplot(data=drawing_data, x=MSP, y=cell, color=color, ax=ax)
    matplotlib.pyplot.text(x=get_middle(drawing_data[MSP]), y=get_middle(drawing_data[cell]), s=f"r={r:.3f}, p={p:.3f}", horizontalalignment="center", verticalalignment="center", fontsize="small", bbox={"alpha": 0.3, "color": "white"})

    matplotlib.pyplot.tight_layout()

    cell_name = cell
    for r in replaced:
        cell_name = cell_name.replace(r, "")

    fig_name = f"Scatter-{stage}-{MSP}-{cell_name}.pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)
    return fig_name


def joint(stage, MSP, cell):
    if stage == "All":
        drawing_data = input_data.copy()
        color = "tab:blue"
    elif stage == "Precancer":
        drawing_data = input_data.loc[~(input_data["Stage"].isin({"Normal", "Primary"}))]
        color = "tab:pink"
    else:
        drawing_data = input_data.loc[(input_data["Stage"] == stage)]
        color = step00.stage_color_code[stage]

    if drawing_data.empty:
        return ""
    elif (numpy.std(drawing_data[MSP]) == 0) and (numpy.std(drawing_data[cell]) == 0):
        return ""

    slope, intercept, r, p, *_ = scipy.stats.linregress(drawing_data[MSP], drawing_data[cell])

    if p >= 0.05:
        return ""

    g = seaborn.jointplot(data=drawing_data, x=MSP, y=cell, color=color, kind="reg", height=24, ratio=5)
    g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})

    cell_name = cell
    for r in replaced:
        cell_name = cell_name.replace(r, "")

    fig_name = f"Joint-{stage}-{MSP}-{cell_name}.pdf"
    g.savefig(fig_name)
    matplotlib.pyplot.close(g.fig)
    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Deconvolution result TSV file (not necessarily TSV)", type=str)
    parser.add_argument("clinical", help="Clinical data w/ Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus <= 1:
        raise ValueError("Number of CPUs must be positive!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    else:
        raise Exception("Something went wrong!!")
    print(clinical_data)

    patients = set(clinical_data.index)
    print(patients)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data = input_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index))), :]
    print(input_data)

    cells = sorted(input_data.columns)

    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    for MSP in tqdm.tqdm(step00.sharing_columns):
        input_data[MSP] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), MSP], list(input_data.index)))
    print(input_data)

    with multiprocessing.Pool(args.cpus) as pool:
        figures = list(pool.starmap(reg, list(itertools.product(step00.long_sample_type_list + ["Precancer", "All"], step00.sharing_columns, cells))))
        figures += list(pool.starmap(joint, list(itertools.product(step00.long_sample_type_list + ["Precancer", "All"], step00.sharing_columns, cells))))

    figures = list(filter(None, figures))

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
