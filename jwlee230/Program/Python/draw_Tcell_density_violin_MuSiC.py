"""
draw_Tcell_density_violin_MuSiC.py: draw T-cell density violin plot for MuSiC
"""
import argparse
import itertools
import matplotlib
import matplotlib.pyplot
import pandas
import scipy.stats
import seaborn
import statannotations.Annotator
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="MuSiC result TSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .tsv!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .pdf!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data = input_data.loc[:, step00.Tcell_list]
    input_data["T cell density"] = list(map(lambda x: sum(input_data.loc[x, step00.Tcell_list]), list(input_data.index)))
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    print(input_data)

    for stage in tqdm.tqdm(set(input_data["Stage"])):
        if len(input_data.loc[(input_data["Stage"] == stage)]) < 3:
            input_data = input_data.loc[~(input_data["Stage"] == stage)]
    print(input_data)

    stage_list = set(input_data["Stage"])
    order = list(filter(lambda x: x in stage_list, step00.long_sample_type_list))
    palette = list(map(lambda x: step00.stage_color_code[x], order))

    stat, p = scipy.stats.kruskal(*[input_data.loc[(input_data["Stage"] == stage), "T cell density"] for stage in order])

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.violinplot(data=input_data, x="Stage", y="T cell density", order=order, palette=palette)
    statannotations.Annotator.Annotator(ax, list(itertools.combinations(order, 2)), data=input_data, x="Stage", y="T cell density", order=order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.title(f"Kruskal-Wallis p={p:.3f}")
    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
