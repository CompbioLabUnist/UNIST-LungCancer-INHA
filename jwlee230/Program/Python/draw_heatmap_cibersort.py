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

    cibersort_data = pandas.read_csv(args.cibersort, sep="\t", index_col="Mixture")
    print(cibersort_data)
    cibersort_data = cibersort_data.reindex(index=sorted(list(cibersort_data.index), key=step00.sorting_by_type)).T
    cibersort_data = cibersort_data.reindex(index=sorted(list(cibersort_data.index), key=lambda x: sum(cibersort_data.loc[x, :]), reverse=True))
    print(cibersort_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(cibersort_data.shape[1], cibersort_data.shape[0]))

    seaborn.heatmap(data=cibersort_data, xticklabels=True, yticklabels=True, cbar=False, square=False, ax=ax, cmap="Reds")

    if args.ADC:
        matplotlib.pyplot.title("ADC")
    elif args.SQC:
        matplotlib.pyplot.title("SQC")
    else:
        raise Exception("Something went wrong!!")
    matplotlib.pyplot.xticks(fontsize="xx-small")
    matplotlib.pyplot.yticks(fontsize="xx-small")
    ax.figure.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
