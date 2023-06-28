"""
draw_clustermap_RSEM_clinical.py: draw clustermap upon RSEM DEG data with clinical information
"""
import argparse
import itertools
import matplotlib
import matplotlib.colors
import matplotlib.patches
import matplotlib.pyplot
import pandas
import seaborn
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="RSEM result TSV file", type=str)
    parser.add_argument("clinical", help="Clinical CSV file", type=str)
    parser.add_argument("DEG", help="DEG result TSV file", type=str, nargs="+")
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--compare", help="Comparison grouping", type=str, default="Recurrence")
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.05)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ADC", help="Draw ADC pathway", action="store_true", default=False)
    group.add_argument("--SQC", help="Draw SQC pathway", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif list(filter(lambda x: not x.endswith(".tsv"), args.DEG)):
        raise ValueError("DEG must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif not (0 < args.p < 1):
        raise ValueError("P-value must be between 0 and 1!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col="gene_name").fillna(0)
    print(input_data)

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    samples = list(input_data.columns)
    if args.SQC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
        samples = list(filter(lambda x: step00.get_patient(x) in histology, samples))
    elif args.ADC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
        samples = list(filter(lambda x: step00.get_patient(x) in histology, samples))
    else:
        raise Exception("Something went wrong!!")
    print(samples)

    input_data = input_data.loc[:, samples]
    print(input_data)

    DEG_list = list()
    for DEG_file in tqdm.tqdm(args.DEG):
        DEG_data = pandas.read_csv(DEG_file, sep="\t", index_col=0)
        DEG_data = DEG_data.loc[(DEG_data["padj"] < args.p)]
        DEG_list.append(set(DEG_data.index))

    input_data = input_data.loc[sorted(set.union(*DEG_list)), sorted(input_data.columns)].dropna()
    print(input_data)
    col_colors = pandas.DataFrame(index=input_data.columns)

    stage_set = set(map(step00.get_long_sample_type, list(input_data.columns)))
    stage_list = list(filter(lambda x: x in stage_set, step00.long_sample_type_list))

    clinical_list = sorted(set(clinical_data[args.compare]))
    clinical_colors = dict(zip(clinical_list, itertools.cycle(matplotlib.colors.XKCD_COLORS)))

    if len(stage_set) > 1:
        col_colors["Stage"] = list(map(step00.get_color_by_type, list(col_colors.index)))
    col_colors[args.compare] = list(map(lambda x: clinical_colors[clinical_data.loc[step00.get_patient(x), args.compare]], list(col_colors.index)))
    print(col_colors)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    g = seaborn.clustermap(data=input_data, figsize=(18, 32), row_cluster=(input_data.shape[0] > 3), col_cluster=(input_data.shape[1] > 3), col_colors=col_colors, xticklabels=False, yticklabels=False, square=False, cmap="coolwarm", z_score=1 if (input_data.shape[0] > 10) else 0, center=0, robust=True, dendrogram_ratio=(0.2, 0.2), cbar_pos=(-0.1, 0.6, 0.05, 0.18))

    matplotlib.pyplot.legend([matplotlib.patches.Patch(facecolor=x) for x in clinical_colors.values()], clinical_colors.keys(), title=args.compare, bbox_to_anchor=(0, 1), bbox_transform=matplotlib.pyplot.gcf().transFigure)

    g.ax_heatmap.set_xlabel(f"{input_data.shape[1]} samples")
    g.ax_heatmap.set_ylabel(f"{input_data.shape[0]} genes")

    g.savefig(args.output)
    g.savefig(args.output.replace(".pdf", ".png"))
