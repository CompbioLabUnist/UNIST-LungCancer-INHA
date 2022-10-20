"""
draw_clustermap_TIMER_MSP.py: draw clustermap plot from TIMER result with clinical information
"""
import argparse
import itertools
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="TIMER result CSV file", type=str)
    parser.add_argument("clinical", help="Clinical information w/ Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

    group_threshold = parser.add_mutually_exclusive_group(required=True)
    group_threshold.add_argument("--median", help="Use median threshold", action="store_true", default=False)
    group_threshold.add_argument("--mean", help="Use mean threshold", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".csv"):
        raise ValueError("Input must end with .CSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    input_data = pandas.read_csv(args.input, index_col=0).T
    cells = list(input_data.columns)
    tools = sorted(set(map(lambda x: x.split("_")[-1], cells)))
    print(input_data)

    patients = set(map(step00.get_patient, list(input_data.index))) & set(clinical_data.index)
    clinical_data = clinical_data.loc[sorted(patients), :]
    input_data = input_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index))), :]
    print(clinical_data)
    print(input_data)

    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    print(input_data)

    MSP_threshold = dict()
    for MSP in tqdm.tqdm(step00.sharing_columns):
        if args.median:
            threshold = numpy.median(clinical_data[MSP])
        elif args.mean:
            threshold = numpy.mean(clinical_data[MSP])

        input_data[MSP] = list(map(lambda x: "Lower" if (clinical_data.loc[step00.get_patient(x), MSP] < threshold) else "Higher", list(input_data.index)))
        MSP_threshold[MSP] = threshold
    print(input_data)

    MSP_order = ["Lower", "Higher"]
    MSP_color_dict = {"Lower": "tab:blue", "Higher": "tab:red"}

    figures = list()
    for MSP, tool in tqdm.tqdm(list(itertools.product(step00.sharing_columns, tools))):
        drawing_data = input_data.loc[:, list(filter(lambda x: x.endswith(tool), cells))]
        drawing_data.columns = list(map(lambda x: "_".join(x.split("_")[:-1]), list(drawing_data.columns)))
        drawing_data = drawing_data.loc[list(filter(lambda x: sum(drawing_data.loc[x, :]), list(drawing_data.index))), list(filter(lambda x: sum(drawing_data.loc[:, x]), list(drawing_data.columns)))]
        for index in list(drawing_data.index):
            drawing_data.loc[index, :] = drawing_data.loc[index, :] / sum(drawing_data.loc[index, :])

        color_data = pandas.DataFrame(index=drawing_data.index)
        color_data["Stage"] = list(map(step00.get_color_by_type, list(color_data.index)))
        color_data["MSP"] = list(map(lambda x: MSP_color_dict[input_data.loc[x, MSP]], list(color_data.index)))

        stage_set = set(map(step00.get_long_sample_type, list(drawing_data.index)))
        stage_list = list(filter(lambda x: x in stage_set, step00.long_sample_type_list))

        drawing_data = drawing_data.T

        g = seaborn.clustermap(data=drawing_data, figsize=(32, 18), row_cluster=True, col_cluster=True, cbar_pos=(-0.05, 0.3, 0.02, 0.5), col_colors=color_data, xticklabels=False, yticklabels=True, square=False, z_score=0, cmap="coolwarm", center=0, robust=True)

        g.ax_heatmap.set_ylabel(f"{drawing_data.shape[0]} cell types")
        g.ax_heatmap.set_xlabel(f"{drawing_data.shape[1]} samples")

        matplotlib.pyplot.legend([matplotlib.patches.Patch(facecolor=step00.stage_color_code[x]) for x in stage_list], stage_list, title="Stage", bbox_to_anchor=(0, 1), bbox_transform=matplotlib.pyplot.gcf().transFigure)

        figures.append(f"{MSP}-{tool}-All.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

    for MSP, stage, tool in tqdm.tqdm(list(itertools.product(step00.sharing_columns, step00.long_sample_type_list, tools))):
        drawing_data = input_data.loc[list(filter(lambda x: step00.get_long_sample_type(x) == stage, list(input_data.index))), list(filter(lambda x: x.endswith(tool), cells))]

        if drawing_data.shape[0] < 2:
            continue

        drawing_data.columns = list(map(lambda x: "_".join(x.split("_")[:-1]), list(drawing_data.columns)))
        drawing_data = drawing_data.loc[list(filter(lambda x: sum(drawing_data.loc[x, :]), list(drawing_data.index))), list(filter(lambda x: sum(drawing_data.loc[:, x]), list(drawing_data.columns)))]
        for index in list(drawing_data.index):
            drawing_data.loc[index, :] = drawing_data.loc[index, :] / sum(drawing_data.loc[index, :])

        color_data = pandas.DataFrame(index=drawing_data.index)
        color_data["MSP"] = list(map(lambda x: MSP_color_dict[input_data.loc[x, MSP]], list(color_data.index)))

        stage_set = set(map(step00.get_long_sample_type, list(drawing_data.index)))
        stage_list = list(filter(lambda x: x in stage_set, step00.long_sample_type_list))

        drawing_data = drawing_data.T

        g = seaborn.clustermap(data=drawing_data, figsize=(32, 18), row_cluster=True, col_cluster=True, cbar_pos=(-0.05, 0.3, 0.02, 0.5), col_colors=color_data, xticklabels=False, yticklabels=True, square=False, z_score=0, cmap="coolwarm", center=0, robust=True)

        g.ax_heatmap.set_ylabel(f"{drawing_data.shape[0]} cell types")
        g.ax_heatmap.set_xlabel(f"{drawing_data.shape[1]} {stage} samples")

        matplotlib.pyplot.legend([matplotlib.patches.Patch(facecolor=MSP_color_dict[x]) for x in MSP_order], MSP_order, title="MSP", bbox_to_anchor=(0, 1), bbox_transform=matplotlib.pyplot.gcf().transFigure)

        figures.append(f"{MSP}-{tool}-{stage}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(figures):
            tar.add(f, arcname=f)
