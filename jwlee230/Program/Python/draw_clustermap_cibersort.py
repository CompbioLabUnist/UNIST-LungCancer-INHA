"""
draw_heatmap_cibersort.py: draw heatmap plot from CibersortX result
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

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ADC", help="Draw ADC pathway", action="store_true", default=False)
    group.add_argument("--SQC", help="Draw SQC pathway", action="store_true", default=False)

    args = parser.parse_args()

    if not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    cibersort_data = pandas.read_csv(args.cibersort, sep="\t", index_col="Mixture").T
    print(cibersort_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    # seaborn.heatmap(data=cibersort_data, xticklabels=True, yticklabels=True, cbar=False, square=False, ax=ax, cmap="coolwarm")
    g = seaborn.clustermap(data=cibersort_data, figsize=(cibersort_data.shape[1], cibersort_data.shape[0]), row_cluster=False, col_cluster=True, cbar_pos=None, col_colors=list(map(step00.get_color_by_type, cibersort_data.columns)), xticklabels=True, yticklabels=True, square=False, cmap="coolwarm")

    g.savefig(args.output)
