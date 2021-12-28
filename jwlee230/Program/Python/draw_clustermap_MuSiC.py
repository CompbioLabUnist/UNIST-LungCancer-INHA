"""
draw_clustermap_MuSiC.py: draw clustermap plot from MuSiC result
"""
import argparse
import matplotlib
import matplotlib.patches
import matplotlib.pyplot
import pandas
import seaborn
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="input result TSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .tsv!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    input_data = pandas.read_csv(args.input, sep="\t")
    input_data = input_data.set_index(list(input_data.columns)[0]).T
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    palette = list(map(lambda x: step00.stage_color_code[step00.get_long_sample_type(x)], list(input_data.columns)))
    stage_set = set(map(step00.get_long_sample_type, list(input_data.columns)))
    stage_list = list(filter(lambda x: x in stage_set, step00.long_sample_type_list))

    g = seaborn.clustermap(data=input_data, figsize=(max(input_data.shape), input_data.shape[0]), row_cluster=True, col_cluster=True, cbar_pos=(-0.04, 0.2, 0.02, 0.6), col_colors=palette, tree_kws={"linewidths": 2.0}, xticklabels=True, yticklabels=True, square=False, cmap="Reds", vmin=0, vmax=1)

    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")

    matplotlib.pyplot.legend([matplotlib.patches.Patch(facecolor=step00.stage_color_code[x]) for x in stage_list], stage_list, title="Stages", bbox_to_anchor=(0, 1), bbox_transform=matplotlib.pyplot.gcf().transFigure)

    g.savefig(args.output)
