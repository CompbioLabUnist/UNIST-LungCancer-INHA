"""
draw_clustermap_BisqueRNA.py: draw clustermap plot from BisqueRNA result
"""
import argparse
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input result TSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .tsv!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    input_data = pandas.read_csv(args.input, sep="\t")
    input_data = input_data.set_index(list(input_data.columns)[0])
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    g = seaborn.clustermap(data=input_data, figsize=(max(input_data.shape), input_data.shape[0]), row_cluster=True, col_cluster=True, cbar_pos=None, col_colors=list(map(step00.get_color_by_type, input_data.columns)), tree_kws={"linewidths": 2.0}, xticklabels=True, yticklabels=True, square=False, cmap="Reds")
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")

    g.savefig(args.output)
