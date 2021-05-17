"""
draw_violin_plots.py: draw violin plots upon DEG
"""
import argparse
import itertools
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import statannot
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("count", help="Count TSV file", type=str)
    parser.add_argument("DEG", help="DEG TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--pvalue", help="P-value threshold", type=float, default=0.05)
    parser.add_argument("--fold", help="Fold change threshold", type=float, default=2)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ADC", help="Draw ADC pathway", action="store_true", default=False)
    group.add_argument("--SQC", help="Draw SQC pathway", action="store_true", default=False)

    args = parser.parse_args()

    if not args.count.endswith(".tsv"):
        raise ValueError("Count must end with .TSV!!")
    elif not args.DEG.endswith(".tsv"):
        raise ValueError("DEG must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    tar_files = list()

    count_data = pandas.read_csv(args.count, sep="\t", index_col="gene_name").T
    count_data["Stage"] = list(map(step00.get_long_sample_type, list(count_data.index)))
    print(count_data)

    DEG_data = pandas.read_csv(args.DEG, sep="\t", header=0, names=["gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"], index_col="gene_name")
    print(DEG_data)

    up_gene = DEG_data.loc[(DEG_data["log2FoldChange"] >= numpy.log2(args.fold)) & (DEG_data["pvalue"] < args.pvalue), :]
    down_gene = DEG_data.loc[(DEG_data["log2FoldChange"] <= -1 * numpy.log2(args.fold)) & (DEG_data["pvalue"] < args.pvalue), :]

    for gene in sorted(list(up_gene.index) + list(down_gene.index)):
        print(gene)

        matplotlib.use("Agg")
        matplotlib.rcParams.update(step00.matplotlib_parameters)
        seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        if args.ADC:
            seaborn.violinplot(data=count_data, x="Stage", y=gene, order=step00.ADC_stage_list)
            statannot.add_stat_annotation(ax, data=count_data, x="Stage", y=gene, order=step00.ADC_stage_list, test="t-test_ind", box_pairs=itertools.combinations(step00.ADC_stage_list, 2), text_format="star", loc="inside", verbose=0)
        elif args.SQC:
            seaborn.violinplot(data=count_data, x="Stage", y=gene, order=step00.SQC_stage_list)
            statannot.add_stat_annotation(ax, data=count_data, x="Stage", y=gene, order=step00.SQC_stage_list, test="t-test_ind", box_pairs=itertools.combinations(step00.SQC_stage_list, 2), text_format="star", loc="inside", verbose=0)
        else:
            raise Exception("Something went wrong!!")

        tar_files.append(gene + ".pdf")
        fig.savefig(tar_files[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for f in tar_files:
            tar.add(f, arcname=f)
