"""
draw_clustermap_RSEM_geneset.py: draw clustermap upon RSEM DEG data with given gene set
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

    parser.add_argument("input", help="RSEM result TSV file", type=str)
    parser.add_argument("clinical", help="Clinical CSV file", type=str)
    parser.add_argument("geneset", help="Gene set TXT file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ADC", help="Draw ADC pathway", action="store_true", default=False)
    group.add_argument("--SQC", help="Draw SQC pathway", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.geneset.endswith(".txt"):
        raise ValueError("GENESET must end with .TXT!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col="gene_name").fillna(0)
    print(input_data)

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    patients = list(input_data.columns)
    if args.SQC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
        patients = list(filter(lambda x: step00.get_patient(x) in histology, patients))
    elif args.ADC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
        patients = list(filter(lambda x: step00.get_patient(x) in histology, patients))
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    input_data = input_data.loc[:, patients]
    print(input_data)

    geneset_data = pandas.read_csv(args.geneset, comment="#")
    DEG_list = sorted(set(geneset_data[list(geneset_data.columns)[0]]) & set(input_data.index))
    print(geneset_data)

    input_data = input_data.loc[DEG_list, :]
    print(input_data)

    palette = list(map(lambda x: step00.stage_color_code[step00.get_long_sample_type(x)], list(input_data.columns)))
    stage_set = set(map(step00.get_long_sample_type, list(input_data.columns)))
    stage_list = list(filter(lambda x: x in stage_set, step00.long_sample_type_list))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    g = seaborn.clustermap(data=input_data, figsize=(18, 32), row_cluster=True, col_cluster=True, col_colors=list(map(step00.get_color_by_type, input_data.columns)), xticklabels=False, yticklabels=False, square=False, cmap="coolwarm", z_score=1, center=0, robust=True)

    matplotlib.pyplot.legend([matplotlib.patches.Patch(facecolor=step00.stage_color_code[x]) for x in stage_list], stage_list, title="Stages", bbox_to_anchor=(0, 1), bbox_transform=matplotlib.pyplot.gcf().transFigure)

    g.ax_heatmap.set_xlabel(f"{input_data.shape[1]} samples")
    g.ax_heatmap.set_ylabel(f"{input_data.shape[0]} genes")

    g.savefig(args.output)
    g.savefig(args.output.replace(".pdf", ".png"))
