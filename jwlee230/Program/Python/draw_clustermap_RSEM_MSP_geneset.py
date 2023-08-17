"""
draw_clustermap_RSEM_MSP_geneset.py: draw clustermap upon RSEM DEG data with MSP with given gene set
"""
import argparse
import matplotlib
import matplotlib.patches
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="RSEM result TSV file", type=str)
    parser.add_argument("clinical", help="Clinical with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("geneset", help="Gene set TXT file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--compare", help="Comparison MSP", type=str, choices=step00.sharing_columns, default=step00.sharing_columns[0])

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ADC", help="Draw ADC pathway", action="store_true", default=False)
    group.add_argument("--SQC", help="Draw SQC pathway", action="store_true", default=False)

    group_threshold = parser.add_mutually_exclusive_group(required=True)
    group_threshold.add_argument("--median", help="Use median threshold", action="store_true", default=False)
    group_threshold.add_argument("--mean", help="Use mean threshold", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.geneset.endswith(".txt"):
        raise ValueError("GENESET must end with .TXT!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col="gene_name").fillna(0)
    print(input_data)

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
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

    geneset_data = pandas.read_csv(args.geneset, comment="#")
    DEG_list = sorted(set(geneset_data[list(geneset_data.columns)[0]]) & set(input_data.index))
    print(geneset_data)

    input_data = input_data.loc[DEG_list, sorted(input_data.columns, key=lambda x: (step00.long_sample_type_list.index(step00.get_long_sample_type(x)), clinical_data.loc[step00.get_patient(x), args.compare]))]
    print(input_data)

    if args.median:
        threshold = numpy.median(clinical_data.loc[list(set(map(step00.get_patient, samples))), args.compare])
    elif args.mean:
        threshold = numpy.mean(clinical_data.loc[list(set(map(step00.get_patient, samples))), args.compare])
    else:
        raise Exception("Something went wrong!!")
    print(args.compare, samples)

    stage_set = set(map(step00.get_long_sample_type, list(input_data.columns)))
    stage_list = list(filter(lambda x: x in stage_set, step00.long_sample_type_list))

    col_colors = pandas.DataFrame(index=input_data.columns)
    if len(stage_set) > 1:
        col_colors["Stage"] = list(map(step00.get_color_by_type, list(col_colors.index)))
    col_colors["MSP"] = list(map(lambda x: "b" if (clinical_data.loc[step00.get_patient(x), args.compare] < threshold) else "r", list(col_colors.index)))
    print(col_colors)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    g = seaborn.clustermap(data=input_data, figsize=(18, 32), row_cluster=(input_data.shape[0] > 1), col_cluster=(len(stage_set) == 1), col_colors=col_colors, xticklabels=False, yticklabels=False, square=False, cmap="coolwarm", z_score=1 if (input_data.shape[0] > 5) else 0, center=0, robust=True, dendrogram_ratio=(0.2 if (input_data.shape[0] > 1) else 0.01, 0.2 if (len(stage_set) == 1) else 0.01), cbar_pos=(-0.1, 0.6, 0.05, 0.18))

    if len(stage_set) > 1:
        matplotlib.pyplot.legend([matplotlib.patches.Patch(facecolor=step00.stage_color_code[x]) for x in stage_list], stage_list, title="Stages", bbox_to_anchor=(0, 1), bbox_transform=matplotlib.pyplot.gcf().transFigure)
    else:
        matplotlib.pyplot.legend([matplotlib.patches.Patch(facecolor=x) for x in ["b", "r"]], ["Lower", "Higher"], title="MSP", bbox_to_anchor=(0, 1), bbox_transform=matplotlib.pyplot.gcf().transFigure)

    g.ax_heatmap.set_xlabel(f"{input_data.shape[1]} samples")
    g.ax_heatmap.set_ylabel(f"{input_data.shape[0]} genes")

    g.savefig(args.output)
    g.savefig(args.output.replace(".pdf", ".png"))
