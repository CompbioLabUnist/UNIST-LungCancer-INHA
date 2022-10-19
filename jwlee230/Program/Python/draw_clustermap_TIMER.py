"""
draw_clustermap_TIMER.py: draw clustermap plot from TIMER result
"""
import argparse
import itertools
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="TIMER result CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".csv"):
        raise ValueError("Input must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    input_data = pandas.read_csv(args.input, index_col=0).T
    cells = list(input_data.columns)
    tools = sorted(set(map(lambda x: x.split("_")[-1], cells)))
    print(input_data)

    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    for stage in set(input_data["Stage"]):
        if len(input_data.loc[(input_data["Stage"] == stage)]) < 3:
            input_data = input_data.loc[~(input_data["Stage"] == stage)]
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    for tool in tqdm.tqdm(tools):
        drawing_data = input_data.loc[:, list(filter(lambda x: x.endswith(tool), cells))]
        drawing_data.columns = list(map(lambda x: "_".join(x.split("_")[:-1]), list(drawing_data.columns)))
        drawing_data = drawing_data.loc[list(filter(lambda x: sum(drawing_data.loc[x, :]), list(drawing_data.index))), list(filter(lambda x: sum(drawing_data.loc[:, x]), list(drawing_data.columns)))]
        for index in list(drawing_data.index):
            drawing_data.loc[index, :] = drawing_data.loc[index, :] / sum(drawing_data.loc[index, :])

        palette = list(map(lambda x: step00.stage_color_code[step00.get_long_sample_type(x)], list(drawing_data.index)))
        stage_set = set(map(step00.get_long_sample_type, list(drawing_data.index)))
        stage_list = list(filter(lambda x: x in stage_set, step00.long_sample_type_list))

        g = seaborn.clustermap(data=drawing_data, figsize=(18, 32), row_cluster=True, col_cluster=True, cbar_pos=(-0.05, 0.3, 0.02, 0.5), row_colors=palette, xticklabels=True, yticklabels=False, square=False, z_score=1, cmap="coolwarm", center=0, robust=True)

        g.ax_heatmap.set_xticklabels(list(drawing_data.columns), fontsize="xx-small")
        g.ax_heatmap.set_xlabel(f"{drawing_data.shape[1]} cell types")
        g.ax_heatmap.set_ylabel(f"{drawing_data.shape[0]} samples")

        matplotlib.pyplot.legend([matplotlib.patches.Patch(facecolor=step00.stage_color_code[x]) for x in stage_list], stage_list, title="Stages", bbox_to_anchor=(0, 1), bbox_transform=matplotlib.pyplot.gcf().transFigure)

        figures.append(f"{tool}-All.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

    for stage, tool in tqdm.tqdm(list(itertools.product(step00.long_sample_type_list, tools))):
        drawing_data = input_data.loc[list(filter(lambda x: step00.get_long_sample_type(x) == stage, list(input_data.index))), list(filter(lambda x: x.endswith(tool), cells))]

        if drawing_data.shape[0] < 2:
            continue

        drawing_data.columns = list(map(lambda x: "_".join(x.split("_")[:-1]), list(drawing_data.columns)))
        drawing_data = drawing_data.loc[list(filter(lambda x: sum(drawing_data.loc[x, :]), list(drawing_data.index))), list(filter(lambda x: sum(drawing_data.loc[:, x]), list(drawing_data.columns)))]
        for index in list(drawing_data.index):
            drawing_data.loc[index, :] = drawing_data.loc[index, :] / sum(drawing_data.loc[index, :])

        g = seaborn.clustermap(data=drawing_data, figsize=(18, 32), row_cluster=True, col_cluster=True, cbar_pos=(-0.05, 0.3, 0.02, 0.5), xticklabels=True, yticklabels=False, square=False, z_score=1, cmap="coolwarm", center=0, robust=True)

        g.ax_heatmap.set_xticklabels(list(drawing_data.columns), fontsize="xx-small")
        g.ax_heatmap.set_xlabel(f"{drawing_data.shape[1]} cell types")
        g.ax_heatmap.set_ylabel(f"{drawing_data.shape[0]} {stage} samples")

        figures.append(f"{tool}-{stage}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(figures):
            tar.add(f, arcname=f)
