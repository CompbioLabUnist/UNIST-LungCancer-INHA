"""
aggregate_PathSeq_clustermap.py: Aggregate PathSeq results as clustermap
"""
import argparse
import itertools
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import pandas
import seaborn
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PathSeq results TSV file", type=str)
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--level", choices=step00.PathSeq_type_list, type=str, required=True)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_sorting = parser.add_mutually_exclusive_group(required=True)
    group_sorting.add_argument("--patient", help="Sorting by patient first", action="store_true", default=False)
    group_sorting.add_argument("--type", help="Sorting by type first", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    output_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    index = list(filter(lambda x: step00.get_patient(x) in patients, list(output_data.index)))

    if args.patient:
        index.sort(key=step00.sorting)
    elif args.type:
        index.sort(key=step00.sorting_by_type)
    else:
        raise Exception("Something went wrong!!")

    output_data = output_data.loc[index, :]
    print(output_data)

    taxa_list = list(output_data.columns)[:-1]
    taxa_coloring = dict(zip(taxa_list, itertools.cycle(matplotlib.colors.XKCD_COLORS)))

    order = list(filter(lambda x: x in set(output_data["Subtype"]), step00.long_sample_type_list))
    print(order)

    row_colors = list(map(lambda x: step00.stage_color_code[step00.get_long_sample_type(x)], list(output_data.index)))
    stage_set = set(map(step00.get_long_sample_type, list(output_data.index)))
    stage_list = list(filter(lambda x: x in stage_set, step00.long_sample_type_list))

    col_colors = list(map(lambda x: taxa_coloring[x], list(output_data.columns)[:-1]))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    try:
        g = seaborn.clustermap(data=output_data.loc[:, taxa_list], figsize=(32, 18), row_cluster=True, col_cluster=True, cbar_pos=(-0.04, 0.2, 0.02, 0.6), row_colors=row_colors, col_colors=col_colors, xticklabels=False, yticklabels=False, square=False, cmap="Reds", vmin=0, vmax=100)
    except RecursionError:
        g = seaborn.clustermap(data=output_data.loc[:, taxa_list], figsize=(32, 18), row_cluster=True, col_cluster=False, cbar_pos=(-0.04, 0.2, 0.02, 0.6), row_colors=row_colors, col_colors=col_colors, xticklabels=False, yticklabels=False, square=False, cmap="Reds", vmin=0, vmax=100, dendrogram_ratio=(0.2, 0.0))

    g.ax_heatmap.set_xlabel(f"{len(taxa_list)} {args.level}")
    g.ax_heatmap.set_ylabel(f"{len(output_data)} samples")

    matplotlib.pyplot.legend([matplotlib.patches.Patch(facecolor=step00.stage_color_code[x]) for x in stage_list], stage_list, title="Stages", bbox_to_anchor=(0, 1), bbox_transform=matplotlib.pyplot.gcf().transFigure)

    g.savefig(args.output)
    g.savefig(args.output.replace(".pdf", ".png"))
