"""
draw_IHC.py: draw IHC images from patient information
"""
import argparse
import itertools
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import statannotations.Annotator
import tqdm
import step00


def change_stage(stage):
    if stage == "Ca":
        return "Primary"
    elif stage == "N":
        return "Normal"
    else:
        return "Precancer"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="IHC result CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

    args = parser.parse_args()

    input_data = pandas.read_csv(args.input)
    print(input_data)

    gene_list = list(input_data.columns)[5:]
    print(len(gene_list), sorted(gene_list))

    input_data["Stage"] = list(map(change_stage, input_data["Dx"]))
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    order = [False, True]
    hue_order = ["Normal", "Precancer", "Primary"]
    palette = {"Normal": "tab:cyan", "Precancer": "tab:green", "Primary": "tab:purple"}

    figures = list()
    for gene in tqdm.tqdm(gene_list):
        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.violinplot(data=input_data, x="Stage", order=hue_order, y=gene, palette=palette, cut=0, linewidth=5, ax=ax)
        if any(numpy.isnan(input_data[gene])):
            statannotations.Annotator.Annotator(ax, list(itertools.combinations(hue_order[1:], r=2)), data=input_data, x="Stage", order=hue_order, y=gene).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()
        else:
            statannotations.Annotator.Annotator(ax, list(itertools.combinations(hue_order, r=2)), data=input_data, x="Stage", order=hue_order, y=gene).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()

        matplotlib.pyplot.ylabel("IHC expression")
        matplotlib.pyplot.title(gene)
        matplotlib.pyplot.tight_layout()

        figures.append(f"{gene}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    for gene in tqdm.tqdm(gene_list):
        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.violinplot(data=input_data, x="Recurrence", order=order, y=gene, hue="Stage", hue_order=hue_order, palette=palette, cut=0, linewidth=5, ax=ax)
        if any(numpy.isnan(input_data[gene])):
            statannotations.Annotator.Annotator(ax, list(filter(lambda x: x[0][1] == x[1][1], itertools.combinations(itertools.product(order, hue_order[1:]), r=2))), data=input_data, x="Recurrence", order=order, y=gene, hue="Stage", hue_order=hue_order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()
        else:
            statannotations.Annotator.Annotator(ax, list(filter(lambda x: x[0][1] == x[1][1], itertools.combinations(itertools.product(order, hue_order), r=2))), data=input_data, x="Recurrence", order=order, y=gene, hue="Stage", hue_order=hue_order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()

        matplotlib.pyplot.ylabel("IHC expression")
        matplotlib.pyplot.title(gene)
        matplotlib.pyplot.legend(loc="lower right")
        matplotlib.pyplot.tight_layout()

        figures.append(f"{gene}_Recurrence.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    for gene in tqdm.tqdm(gene_list):
        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.violinplot(data=input_data, x="Recurrence", order=order, y=gene, hue="Stage", hue_order=hue_order, palette=palette, cut=0, linewidth=5, ax=ax)
        if any(numpy.isnan(input_data[gene])):
            statannotations.Annotator.Annotator(ax, list(filter(lambda x: x[0][0] == x[1][0], itertools.combinations(itertools.product(order, hue_order[1:]), r=2))), data=input_data, x="Recurrence", order=order, y=gene, hue="Stage", hue_order=hue_order,).configure(test="Mann-Whitney", text_format="star", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()
        else:
            statannotations.Annotator.Annotator(ax, list(filter(lambda x: x[0][0] == x[1][0], itertools.combinations(itertools.product(order, hue_order), r=2))), data=input_data, x="Recurrence", order=order, y=gene, hue="Stage", hue_order=hue_order,).configure(test="Mann-Whitney", text_format="star", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()

        matplotlib.pyplot.ylabel("IHC expression")
        matplotlib.pyplot.title(gene)
        matplotlib.pyplot.legend(loc="lower right")
        matplotlib.pyplot.tight_layout()

        figures.append(f"{gene}_Stage.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(figures):
            tar.add(f, arcname=f)
