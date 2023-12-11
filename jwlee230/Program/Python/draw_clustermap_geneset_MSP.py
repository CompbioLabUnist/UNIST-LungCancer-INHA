"""
draw_clustermap_geneset_MSP.py: draw clustermap with gene set on MSP
"""
import argparse
import itertools
import tarfile
import typing
import matplotlib
import matplotlib.patches
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="RSEM result TSV file", type=str)
    parser.add_argument("clinical", help="Clinical with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("geneset", help="Gene set TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ADC", help="Draw ADC pathway", action="store_true", default=False)
    group.add_argument("--SQC", help="Draw SQC pathway", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.geneset.endswith(".tsv"):
        raise ValueError("Gene-set must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col="gene_name").fillna(0)
    print(input_data)

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients: typing.Set[str] = set(clinical_data.index)
    print(patients)

    geneset_data = pandas.read_csv(args.geneset, sep="\t", index_col=0, header=0, names=["Index", "Contents"])
    print(geneset_data)

    notable_genes: typing.List[str] = sorted(list(filter(None, geneset_data.loc["GENE_SYMBOLS", "Contents"].split(","))))
    print(len(notable_genes), notable_genes)

    input_data = input_data.loc[list(filter(lambda x: x in set(input_data.index), notable_genes)), list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.columns)))]
    input_data = input_data.loc[list(filter(lambda x: numpy.std(input_data.loc[x, :]) > 0, list(input_data.index))), :]
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    for stage, MSP in tqdm.tqdm(list(itertools.product(step00.long_sample_type_list, step00.sharing_columns))):
        drawing_data = input_data.loc[:, sorted(list(filter(lambda x: step00.get_long_sample_type(x) == stage, list(input_data.columns))), key=lambda x: clinical_data.loc[step00.get_patient(x), MSP], reverse=True)]
        drawing_data = drawing_data.loc[list(filter(lambda x: numpy.std(drawing_data.loc[x, :]) > 0, list(drawing_data.index))), :]

        if drawing_data.empty:
            continue

        MSP_list = list(map(lambda x: clinical_data.loc[step00.get_patient(x), MSP], list(drawing_data.columns)))
        col_colors = list(map(lambda x: ((max(MSP_list) - x) / max(MSP_list), (max(MSP_list) - x) / max(MSP_list), (max(MSP_list) - x) / max(MSP_list)), MSP_list))

        g = seaborn.clustermap(data=drawing_data, figsize=(18, 32), row_cluster=True, col_cluster=False, col_colors=col_colors, xticklabels=False, yticklabels=False, square=False, cmap="coolwarm", z_score=0, center=0, robust=True)

        g.ax_heatmap.set_xlabel(f"{drawing_data.shape[1]} {stage} samples")
        g.ax_heatmap.set_ylabel(f"{drawing_data.shape[0]} genes")

        figures.append(f"{stage}-{MSP}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

    for MSP in tqdm.tqdm(step00.sharing_columns):
        input_data = input_data.loc[:, sorted(list(input_data.columns), key=lambda x: clinical_data.loc[step00.get_patient(x), MSP], reverse=True)]

        MSP_list = list(map(lambda x: clinical_data.loc[step00.get_patient(x), MSP], list(input_data.columns)))
        col_colors = list(map(lambda x: ((max(MSP_list) - x) / max(MSP_list), (max(MSP_list) - x) / max(MSP_list), (max(MSP_list) - x) / max(MSP_list)), MSP_list))

        g = seaborn.clustermap(data=input_data, figsize=(18, 32), row_cluster=True, col_cluster=False, col_colors=col_colors, xticklabels=False, yticklabels=False, square=False, cmap="coolwarm", z_score=0, center=0, robust=True)

        g.ax_heatmap.set_xlabel(f"{input_data.shape[1]} samples")
        g.ax_heatmap.set_ylabel(f"{input_data.shape[0]} genes")

        figures.append(f"All-{MSP}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(figures):
            tar.add(f, arcname=f)
