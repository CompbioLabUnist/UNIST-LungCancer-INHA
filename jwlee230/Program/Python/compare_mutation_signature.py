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


def run(stage: str, signature: str, column: str) -> str:
    tmp_data = signature_data.loc[(signature_data["Stage"] == stage), [column, signature]]
    if tmp_data.shape[0] < 3:
        return ""

    r, p = scipy.stats.pearsonr(tmp_data[signature], tmp_data[column])

    g = seaborn.jointplot(data=tmp_data, x=signature, y=column, kind="reg", height=24, ratio=6, color=step00.stage_color_code[stage])
    g.fig.text(0.5, 0.75, "r={0:.3f}, p={1:.3f}".format(r, p), color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")
    g.set_axis_labels("{0} proportion in {1}".format(signature, stage), column)

    fig_name = f"{stage}-{column}-{signature}.pdf"
    g.savefig(fig_name)
    matplotlib.pyplot.close(g.fig)

    return fig_name


def run_all(signature: str, column: str) -> str:
    tmp_data = signature_data.loc[:, [column, signature]]

    if tmp_data.shape[0] < 3:
        return ""

    r, p = scipy.stats.pearsonr(tmp_data[signature], tmp_data[column])

    g = seaborn.jointplot(data=tmp_data, x=signature, y=column, kind="reg", height=24, ratio=6)
    g.fig.text(0.5, 0.75, "r={0:.3f}, p={1:.3f}".format(r, p), color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")
    g.plot_marginals(seaborn.histplot, kde=True, stat="probability", multiple="stack")
    g.set_axis_labels("{0} proportion".format(signature), column)

    fig_name = f"All-{column}-{signature}.pdf"
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
        figures += pool.starmap(run, itertools.product(step00.long_sample_type_list, signature_list, step00.sharing_columns))
        figures += pool.starmap(run_all, itertools.product(signature_list, step00.sharing_columns))

    figures = list(filter(None, figures))

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
