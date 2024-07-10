"""
compare_mutation_signature.py: compare mutation shared proportion vs. mutational signature from SigProfiler
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

signature_data = pandas.DataFrame()


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


def get_middle(values):
    return (min(values) + max(values)) / 2


def joint(stage: str, signature: str, column: str) -> str:
    tmp_data = signature_data.loc[(signature_data["Stage"] == stage), [column, signature]]
    if tmp_data.shape[0] < 3:
        return ""

    r, p = scipy.stats.pearsonr(tmp_data[column], tmp_data[signature])

    g = seaborn.jointplot(data=tmp_data, x=column, y=signature, kind="reg", height=18, ratio=6, color=step00.stage_color_code[stage])
    g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
    g.set_axis_labels(column, f"{signature} in {stage}")

    fig_name = f"Joint-{stage}-{column}-{signature}.pdf"
    g.savefig(fig_name)
    matplotlib.pyplot.close(g.fig)

    return fig_name


def reg(stage: str, signature: str, column: str) -> str:
    tmp_data = signature_data.loc[(signature_data["Stage"] == stage), [column, signature]]
    if tmp_data.shape[0] < 3:
        return ""

    r, p = scipy.stats.pearsonr(tmp_data[column], tmp_data[signature])

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    seaborn.regplot(data=tmp_data, x=column, y=signature, fit_reg=True, scatter=True, color=step00.stage_color_code[stage], ax=ax)

    matplotlib.pyplot.ylabel(f"{signature} in {stage}")
    matplotlib.pyplot.text(get_middle(tmp_data[column]), get_middle(tmp_data[signature]), f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
    matplotlib.pyplot.tight_layout()

    fig_name = f"Scatter-{stage}-{column}-{signature}.pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


def joint_precancer(signature: str, column: str) -> str:
    tmp_data = signature_data.loc[~(signature_data["Stage"].isin({"Normal", "Primary"})), [column, signature]]
    if tmp_data.shape[0] < 3:
        return ""

    r, p = scipy.stats.pearsonr(tmp_data[column], tmp_data[signature])

    g = seaborn.jointplot(data=tmp_data, x=column, y=signature, kind="reg", height=18, ratio=6, color=step00.precancer_color_code["Precancer"])
    g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
    g.plot_marginals(seaborn.histplot, kde=True, stat="probability", multiple="stack")
    g.set_axis_labels(column, f"{signature} in Precancer")

    fig_name = f"Joint-Precancer-{column}-{signature}.pdf"
    g.savefig(fig_name)
    matplotlib.pyplot.close(g.fig)

    return fig_name


def reg_precancer(signature: str, column: str) -> str:
    tmp_data = signature_data.loc[~(signature_data["Stage"].isin({"Normal", "Primary"})), [column, signature]]
    if tmp_data.shape[0] < 3:
        return ""

    r, p = scipy.stats.pearsonr(tmp_data[column], tmp_data[signature])

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    seaborn.regplot(data=tmp_data, x=column, y=signature, fit_reg=True, scatter=True, color=step00.precancer_color_code["Precancer"], ax=ax)

    matplotlib.pyplot.ylabel(f"{signature} in Precancer")
    matplotlib.pyplot.text(get_middle(tmp_data[column]), get_middle(tmp_data[signature]), f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
    matplotlib.pyplot.tight_layout()

    fig_name = f"Scatter-Precancer-{column}-{signature}.pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


def joint_all(signature: str, column: str) -> str:
    tmp_data = signature_data.loc[:, [column, signature]]
    if tmp_data.shape[0] < 3:
        return ""

    r, p = scipy.stats.pearsonr(tmp_data[column], tmp_data[signature])

    g = seaborn.jointplot(data=tmp_data, x=column, y=signature, kind="reg", height=18, ratio=6, color="tab:blue")
    g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
    g.plot_marginals(seaborn.histplot, kde=True, stat="probability", multiple="stack")
    g.set_axis_labels(column, f"{signature}")

    fig_name = f"Joint-All-{column}-{signature}.pdf"
    g.savefig(fig_name)
    matplotlib.pyplot.close(g.fig)

    return fig_name


def reg_all(signature: str, column: str) -> str:
    tmp_data = signature_data.loc[:, [column, signature]]
    if tmp_data.shape[0] < 3:
        return ""

    r, p = scipy.stats.pearsonr(tmp_data[column], tmp_data[signature])

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    seaborn.regplot(data=tmp_data, x=column, y=signature, fit_reg=True, scatter=True, color="tab:blue", ax=ax)

    matplotlib.pyplot.ylabel(f"{signature}")
    matplotlib.pyplot.text(get_middle(tmp_data[column]), get_middle(tmp_data[signature]), f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
    matplotlib.pyplot.tight_layout()

    fig_name = f"Scatter-All-{column}-{signature}.pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


def lm(signature: str, column: str) -> str:
    tmp_data = signature_data.loc[:, [column, signature, "Stage"]]
    tmp_data["Stage"] = list(map(lambda x: x if (x in {"Normal", "Primary"}) else "Precancer", tmp_data["Stage"]))
    if tmp_data.shape[0] < 3:
        return ""

    stages = ["Precancer", "Primary"]
    text = ""
    for stage in stages:
        if tmp_data.loc[(tmp_data["Stage"] == stage), column].shape[0] < 3:
            continue

        r, p = scipy.stats.pearsonr(tmp_data.loc[(tmp_data["Stage"] == stage), column], tmp_data.loc[(tmp_data["Stage"] == stage), signature])
        text += f"{stage}: r={r:.3f}, p={p:.3f}\n"
    text = text.strip()

    g = seaborn.lmplot(data=tmp_data, x=column, y=signature, hue="Stage", hue_order=stages, palette=step00.precancer_color_code, height=18, aspect=1, legend=True, legend_out=False, scatter=True, fit_reg=True)
    g.fig.text(0.5, 0.5, text, color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
    g.set_axis_labels(column, f"{signature}")

    fig_name = f"lm-All-{column}-{signature}.pdf"
    g.savefig(fig_name)
    matplotlib.pyplot.close(g.fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("signature", help="Mutation signature TSV file (not necessarily TSV)", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_relative = parser.add_mutually_exclusive_group(required=True)
    group_relative.add_argument("--absolute", help="Use Absolute(Count) value", action="store_true", default=False)
    group_relative.add_argument("--relative", help="Use Relative(Proportion) value", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    clinical_data: pandas.DataFrame = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(patients)

    signature_data = pandas.read_csv(args.signature, sep="\t", index_col="Samples")
    signature_list = list(signature_data.columns)
    if args.relative:
        for index in tqdm.tqdm(list(signature_data.index)):
            signature_data.loc[index, :] = signature_data.loc[index, :] / sum(signature_data.loc[index, :])
    signature_data["Patient"] = list(map(step00.get_patient, list(signature_data.index)))
    signature_data["Stage"] = list(map(step00.get_long_sample_type, list(signature_data.index)))
    print(signature_data)

    patients &= set(signature_data["Patient"])

    signature_data = signature_data.loc[signature_data["Patient"].isin(patients)]
    for column in tqdm.tqdm(step00.sharing_columns):
        signature_data[column] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), column], list(signature_data.index)))
    print(signature_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()

    with multiprocessing.Pool(args.cpus) as pool:
        figures += pool.starmap(joint, itertools.product(step00.long_sample_type_list, signature_list, step00.sharing_columns))
        figures += pool.starmap(reg, itertools.product(step00.long_sample_type_list, signature_list, step00.sharing_columns))
        figures += pool.starmap(joint_precancer, itertools.product(signature_list, step00.sharing_columns))
        figures += pool.starmap(reg_precancer, itertools.product(signature_list, step00.sharing_columns))
        figures += pool.starmap(joint_all, itertools.product(signature_list, step00.sharing_columns))
        figures += pool.starmap(reg_all, itertools.product(signature_list, step00.sharing_columns))
        figures += pool.starmap(lm, itertools.product(signature_list, step00.sharing_columns))

    figures = list(filter(None, figures))

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
