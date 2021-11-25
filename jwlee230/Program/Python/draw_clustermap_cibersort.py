"""
draw_clustermap_cibersort.py: draw clustermap plot from CibersortX result
"""
import argparse
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("cibersort", help="CIBERSORT result TSV file (not necessarily TSV)", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    args = parser.parse_args()

    if not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    cibersort_data = pandas.read_csv(args.cibersort, sep="\t", index_col="Mixture").T
    print(cibersort_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    g = seaborn.clustermap(data=cibersort_data, figsize=(max(cibersort_data.shape), cibersort_data.shape[0]), row_cluster=True, col_cluster=True, cbar_pos=None, col_colors=list(map(step00.get_color_by_type, cibersort_data.columns)), tree_kws={"linewidths": 2.0}, xticklabels=True, yticklabels=True, square=False, cmap="Reds", vmin=0, vmax=1)
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")

    g.savefig(args.output)
