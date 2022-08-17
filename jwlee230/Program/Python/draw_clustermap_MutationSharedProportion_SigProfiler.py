"""
draw_clustermap_MutationSharedProportion_SigProfiler.py: draw Cluster Map with Mutation Shared Proportion vs. SigProfiler
"""
import argparse
import itertools
import tarfile
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("signature", help="Signature TXT (not necessarily TSV) file", type=str)
    parser.add_argument("clinical", help="Clinical data w/ Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_threshold = parser.add_mutually_exclusive_group(required=True)
    group_threshold.add_argument("--median", help="Use Median threshold", action="store_true", default=False)
    group_threshold.add_argument("--mean", help="Use Mean threshold", action="store_true", default=False)

    group_relative = parser.add_mutually_exclusive_group(required=True)
    group_relative.add_argument("--absolute", help="Use Absolute(Count) value", action="store_true", default=False)
    group_relative.add_argument("--relative", help="Use Relative(Proportion) value", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical data must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    print(clinical_data)

    input_data = pandas.read_csv(args.signature, sep="\t", index_col=0)
    if args.relative:
        for index in tqdm.tqdm(list(input_data.index)):
            input_data.loc[index, :] = input_data.loc[index, :] / sum(input_data.loc[index, :])
    print(input_data)

    signatures = list(input_data.columns)
    input_data["Patient"] = list(map(step00.get_patient, list(input_data.index)))
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    print(input_data)

    input_data = input_data.loc[(input_data["Patient"].isin(set(clinical_data.index)))]
    print(input_data)

    stage_list = list(filter(lambda x: x in set(input_data["Stage"]), step00.long_sample_type_list))
    print(stage_list)

    patient_colors = dict(zip(sorted(set(input_data["Patient"])), itertools.cycle(matplotlib.colors.XKCD_COLORS)))

    figures = list()
    for stage, MSP in tqdm.tqdm(list(itertools.product(stage_list, step00.sharing_columns))):
        drawing_data = input_data.loc[(input_data["Stage"] == stage), signatures]

        if args.median:
            threshold = numpy.median(clinical_data[MSP])
        elif args.mean:
            threshold = numpy.mean(clinical_data[MSP])
        else:
            raise Exception("Something went wrong!!")

        colors = pandas.DataFrame(index=drawing_data.index)
        colors["MSP"] = list(map(lambda x: "tab:blue" if (clinical_data.loc[step00.get_patient(x), MSP] < threshold) else "tab:red", list(drawing_data.index)))

        if args.absolute:
            g = seaborn.clustermap(data=drawing_data, figsize=(18, 32), row_cluster=True, col_cluster=True, cbar_pos=(-0.04, 0.2, 0.02, 0.6), row_colors=colors, tree_kws={"linewidths": 2.0}, xticklabels=True, yticklabels=False, square=False, cmap="Reds")
        elif args.relative:
            g = seaborn.clustermap(data=drawing_data, figsize=(18, 32), row_cluster=True, col_cluster=True, cbar_pos=(-0.04, 0.2, 0.02, 0.6), row_colors=colors, tree_kws={"linewidths": 2.0}, xticklabels=True, yticklabels=False, square=False, cmap="Reds", vmin=0, vmax=1)
        else:
            raise Exception("Something went wrong!!")

        g.ax_heatmap.set_xlabel("")
        g.ax_heatmap.set_ylabel(f"{stage}")

        figures.append(f"{stage}_{MSP}.pdf")
        g.savefig(figures[-1])

    for MSP in tqdm.tqdm(step00.sharing_columns):
        colors = pandas.DataFrame(index=input_data.index)
        colors["Patient"] = list(map(lambda x: patient_colors[step00.get_patient(x)], list(input_data.index)))
        colors["MSP"] = list(map(lambda x: "tab:blue" if (clinical_data.loc[step00.get_patient(x), MSP] < threshold) else "tab:red", list(input_data.index)))

        if args.absolute:
            g = seaborn.clustermap(data=input_data[signatures], figsize=(18, 32), row_cluster=True, col_cluster=True, cbar_pos=(-0.04, 0.2, 0.02, 0.6), row_colors=colors, tree_kws={"linewidths": 2.0}, xticklabels=True, yticklabels=False, square=False, cmap="Reds")
        elif args.relative:
            g = seaborn.clustermap(data=input_data[signatures], figsize=(18, 32), row_cluster=True, col_cluster=True, cbar_pos=(-0.04, 0.2, 0.02, 0.6), row_colors=colors, tree_kws={"linewidths": 2.0}, xticklabels=True, yticklabels=False, square=False, cmap="Reds", vmin=0, vmax=1)
        else:
            raise Exception("Something went wrong!!")

        g.ax_heatmap.set_xlabel("")
        g.ax_heatmap.set_ylabel("All")

        figures.append(f"All_{MSP}.pdf")
        g.savefig(figures[-1])

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
