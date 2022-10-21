"""
draw_violin_TIMER.py: Draw violin plot from TIMER
"""
import argparse
import itertools
import tarfile
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

    parser.add_argument("input", help="TIMER result CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".csv"):
        raise ValueError("Input must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = pandas.read_csv(args.input, index_col=0).T
    cells = list(input_data.columns)
    print(input_data)

    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    for stage in set(input_data["Stage"]):
        if len(input_data.loc[(input_data["Stage"] == stage)]) < 1:
            input_data = input_data.loc[~(input_data["Stage"] == stage)]
    print(input_data)

    tmp = set(input_data["Stage"])
    order = list(filter(lambda x: x in tmp, step00.long_sample_type_list))
    palette = list(map(lambda x: step00.stage_color_code[x], order))

    figures = list()
    for cell in tqdm.tqdm(cells):
        title, tool = cell.split("_")

        try:
            stat, p = scipy.stats.kruskal(*[input_data.loc[(input_data["Stage"] == stage), cell] for stage in order])
        except ValueError:
            _, p = 0.0, 1.0

        compare_list = list()
        for a, b in itertools.combinations(order, r=2):
            try:
                stat, pvalue = scipy.stats.mannwhitneyu(input_data.loc[(input_data["Stage"] == a), cell], input_data.loc[(input_data["Stage"] == b), cell])
            except ValueError:
                continue

            if pvalue < 0.05:
                compare_list.append((a, b))

        if (p > 0.05) and (not compare_list):
            continue

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        seaborn.violinplot(data=input_data, x="Stage", y=cell, order=order, palette=palette, cut=1, linewidth=5, ax=ax)
        if compare_list:
            statannotations.Annotator.Annotator(ax, compare_list, data=input_data, x="Stage", y=cell, order=order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

        matplotlib.pyplot.title(f"Kruskal-Wallis p={p:.3f}")
        matplotlib.pyplot.ylabel(f"{title} from {tool}")
        matplotlib.pyplot.tight_layout()

        figures.append(f"{tool}-{title}.pdf".replace("/", "_"))
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(figures):
            tar.add(f, arcname=f)
