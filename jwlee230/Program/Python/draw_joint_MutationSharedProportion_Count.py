"""
draw_joint_MutationSharedProportion_Count.py: draw joint plot with Mutation Shared Proportion and Mutation Count
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


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


def get_middle(array: list) -> float:
    return (min(array) + max(array)) / 2.0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input MAF files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(patients)

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))
    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
        mutect_data["Tumor_Sample_Barcode"] = pool.map(step00.get_id, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Patient"] = pool.map(step00.get_patient, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Stage"] = pool.map(step00.get_long_sample_type, mutect_data["Tumor_Sample_Barcode"])
    mutect_data = mutect_data.loc[(mutect_data[step00.nonsynonymous_column].isin(step00.nonsynonymous_mutations))]
    print(mutect_data)

    patients &= set(mutect_data["Patient"])

    stage_list = list(filter(lambda x: x in set(mutect_data["Stage"]), step00.long_sample_type_list))
    print(stage_list)

    for stage in tqdm.tqdm(stage_list):
        clinical_data[stage] = None

    for input_data in tqdm.tqdm(args.input):
        clinical_data.loc[step00.get_patient(input_data), step00.get_long_sample_type(input_data)] = mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == step00.get_id(input_data))].shape[0]
    print(clinical_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    for stage, column in tqdm.tqdm(list(itertools.product(stage_list, step00.sharing_columns))):
        drawing_data = clinical_data.dropna(axis="index", subset=[column, stage]).copy()
        drawing_data[stage] = list(map(int, drawing_data[stage]))

        r, p = scipy.stats.pearsonr(drawing_data[column], drawing_data[stage])

        g = seaborn.jointplot(data=drawing_data, x=column, y=stage, kind="reg", height=18, dropna=True, color=step00.stage_color_code[stage])
        g.fig.text(0.5, 0.75, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        g.set_axis_labels(column, f"Mutation Count ({stage})")

        figures.append(f"Joint_{stage}_{column}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

    for stage, column in tqdm.tqdm(list(itertools.product(stage_list, step00.sharing_columns))):
        drawing_data = clinical_data.dropna(axis="index", subset=[column, stage]).copy()
        drawing_data[stage] = list(map(int, drawing_data[stage]))

        r, p = scipy.stats.pearsonr(drawing_data[column], drawing_data[stage])

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.regplot(data=drawing_data, x=column, y=stage, scatter=True, fit_reg=True, color=step00.stage_color_code[stage], ax=ax)
        matplotlib.pyplot.text(get_middle(drawing_data[column]), get_middle(drawing_data[stage]), f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.xlabel(column)
        matplotlib.pyplot.ylabel(f"Mutation Count ({stage})")
        matplotlib.pyplot.tight_layout()

        figures.append(f"Scatter_{stage}_{column}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
