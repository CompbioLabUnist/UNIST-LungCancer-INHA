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

    parser.add_argument("count", help="Count TSV file", type=str)
    parser.add_argument("cibersort", help="CIBERSORT result TSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ADC", help="Draw ADC pathway", action="store_true", default=False)
    group.add_argument("--SQC", help="Draw SQC pathway", action="store_true", default=False)

    args = parser.parse_args()

    if not args.count.endswith(".tsv"):
        raise ValueError("Count must end with .TSV!!")
    elif not args.cibersort.endswith(".tsv"):
        raise ValueError("CIBERSORT must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    count_data = pandas.read_csv(args.count, sep="\t", index_col="gene_name")
    print(count_data)

    cibersort_data = pandas.read_csv(args.cibersort, sep="\t", index_col="Column")
    cibersort_data["Patient_ID"] = list(count_data.columns)
    cibersort_data.set_index("Patient_ID", inplace=True, verify_integrity=True)
    cibersort_data.drop(columns=list(cibersort_data.columns)[-4:], inplace=True)
    cibersort_data = cibersort_data.reindex(index=sorted(list(cibersort_data.index), key=step00.sorting_by_type)).T
    cibersort_data = cibersort_data.reindex(index=sorted(list(cibersort_data.index)))
    print(cibersort_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(cibersort_data.shape[1], cibersort_data.shape[0]))

    seaborn.heatmap(data=cibersort_data, xticklabels=True, yticklabels=True, cbar=False, square=False, ax=ax, cmap="coolwarm")

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
