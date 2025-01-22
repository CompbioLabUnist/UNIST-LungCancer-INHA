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
import scipy.stats
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


def get_middle(values):
    return (max(values) + min(values)) / 2


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="IHC result CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

    args = parser.parse_args()

    input_data = pandas.read_csv(args.input, comment="#")
    print(input_data)

    gene_list = list(input_data.columns)[5:-2]
    print(len(gene_list), sorted(gene_list))

    input_data["Stage"] = list(map(change_stage, input_data["Dx"]))
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    order = [False, True]
    hue_order = ["Precancer", "Primary"]
    short_dict = {"Precancer": "PRE", "Primary": "PRI"}
    palette = {"Precancer": "tab:green", "Primary": "tab:purple"}

    figures = list()
    for gene, stage in tqdm.tqdm(list(itertools.product(gene_list, hue_order))):
        gene_data = input_data.loc[(input_data["Stage"] == stage)].dropna(subset=gene)
        value_list = sorted(set(gene_data[gene]))

        if len(value_list) < 2:
            continue

        count_data = pandas.DataFrame(data=numpy.zeros((len(order), len(value_list))), index=order, columns=value_list, dtype=int)
        for index, row in gene_data.iterrows():
            count_data.loc[row["Recurrence"], row[gene]] += 1

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.heatmap(count_data, cmap="Reds", annot=True, fmt="d", cbar=False, square=False, ax=ax)

        matplotlib.pyplot.title(f"{gene} in {stage}")
        matplotlib.pyplot.xlabel(f"Strength of {gene}")
        matplotlib.pyplot.ylabel("Recurrence")
        matplotlib.pyplot.tight_layout()

        figures.append(f"{gene}_{stage}_Count.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

        for threshold in tqdm.tqdm(value_list[:-1], leave=False):
            merged_count_data = pandas.DataFrame([[0, 0], [0, 0]], index=order, columns=["Low", "High"], dtype=int)

            for index, column in itertools.product(order, value_list):
                if column <= threshold:
                    merged_count_data.loc[index, "Low"] += count_data.loc[index, column]
                elif column > threshold:
                    merged_count_data.loc[index, "High"] += count_data.loc[index, column]
                else:
                    raise ValueError(f"Invalid value: {column}!!")

            p_value = scipy.stats.barnard_exact(merged_count_data.to_numpy()).pvalue

            fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

            seaborn.heatmap(merged_count_data, cmap="Reds", annot=True, fmt="d", cbar=False, square=False, ax=ax)

            matplotlib.pyplot.title(f"{gene} in {short_dict[stage]} (p={p_value:.3f})")
            matplotlib.pyplot.xlabel(f"Strength of {gene} (threshold={int(threshold)})")
            matplotlib.pyplot.ylabel("Recurrence")
            matplotlib.pyplot.tight_layout()

            figures.append(f"{gene}_{stage}_{int(threshold)}.pdf")
            fig.savefig(figures[-1])
            matplotlib.pyplot.close(fig)

    for gene in tqdm.tqdm(gene_list):
        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.swarmplot(data=input_data, x="Recurrence", order=order, y=gene, hue="Stage", hue_order=hue_order, palette=palette, dodge=True, size=40, ax=ax)
        if any(numpy.isnan(input_data[gene])):
            statannotations.Annotator.Annotator(ax, list(filter(lambda x: x[0][1] == x[1][1], itertools.combinations(itertools.product(order, hue_order[1:]), r=2))), data=input_data, x="Recurrence", order=order, y=gene, hue="Stage", hue_order=hue_order).configure(test="Mann-Whitney", text_format="star", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()
        else:
            statannotations.Annotator.Annotator(ax, list(filter(lambda x: x[0][1] == x[1][1], itertools.combinations(itertools.product(order, hue_order), r=2))), data=input_data, x="Recurrence", order=order, y=gene, hue="Stage", hue_order=hue_order).configure(test="Mann-Whitney", text_format="star", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()

        matplotlib.pyplot.ylabel(f"{gene} expression")
        matplotlib.pyplot.title(gene)
        matplotlib.pyplot.legend(loc="lower right")
        matplotlib.pyplot.tight_layout()

        figures.append(f"{gene}_Recurrence.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    for gene in tqdm.tqdm(gene_list):
        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.swarmplot(data=input_data, x="Recurrence", order=order, y=gene, hue="Stage", hue_order=hue_order, palette=palette, dodge=True, size=40, ax=ax)
        if any(numpy.isnan(input_data[gene])):
            statannotations.Annotator.Annotator(ax, list(filter(lambda x: x[0][0] == x[1][0], itertools.combinations(itertools.product(order, hue_order[1:]), r=2))), data=input_data, x="Recurrence", order=order, y=gene, hue="Stage", hue_order=hue_order,).configure(test="Wilcoxon", text_format="star", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()
        else:
            statannotations.Annotator.Annotator(ax, list(filter(lambda x: x[0][0] == x[1][0], itertools.combinations(itertools.product(order, hue_order), r=2))), data=input_data, x="Recurrence", order=order, y=gene, hue="Stage", hue_order=hue_order,).configure(test="Wilcoxon", text_format="star", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()

        matplotlib.pyplot.ylabel(f"{gene} expression")
        matplotlib.pyplot.title(gene)
        matplotlib.pyplot.legend(loc="lower right")
        matplotlib.pyplot.tight_layout()

        figures.append(f"{gene}_Stage.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    for gene, stage, survival in tqdm.tqdm(list(itertools.product(gene_list, hue_order, ["RFS", "OS"]))):
        gene_data = input_data.loc[(input_data["Stage"] == stage)].dropna(subset=gene)

        if gene_data.empty:
            continue

        r_value, p_value = scipy.stats.spearmanr(gene_data[survival], gene_data[gene])

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.regplot(data=gene_data, x=survival, y=gene, scatter=True, fit_reg=True, truncate=False, color=palette[stage], ax=ax)
        matplotlib.pyplot.text(get_middle(gene_data["RFS"]), get_middle(gene_data[gene]), f"r={r_value:.3f}, p={p_value:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center")

        matplotlib.pyplot.title(f"{gene} in {stage}")
        matplotlib.pyplot.xlabel(f"{survival} (days)")
        matplotlib.pyplot.ylabel(f"{gene} expression")
        matplotlib.pyplot.tight_layout()

        figures.append(f"{gene}_{stage}_{survival}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(figures):
            tar.add(f, arcname=f)
