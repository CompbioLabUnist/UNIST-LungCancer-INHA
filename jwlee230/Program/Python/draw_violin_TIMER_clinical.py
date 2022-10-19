"""
draw_violin_TIMER_clinical.py: Draw violin plot from TIMER with clinical information
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
    parser.add_argument("clinical", help="Clinical data data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--column", help="Clinical information column", nargs="+", default=["Recurrence", "NO", "YES"])

    args = parser.parse_args()

    if not args.input.endswith(".csv"):
        raise ValueError("Input must end with .CSV!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    histology = set(clinical_data.index)
    print(histology)

    input_data = pandas.read_csv(args.input, index_col=0).T
    cells = list(input_data.columns)
    print(input_data)

    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    input_data[args.column[0]] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), args.column[0]], list(input_data.index)))
    for stage in set(input_data["Stage"]):
        if any([(len(input_data.loc[(input_data["Stage"] == stage) & (input_data[args.column[0]] == clinical)]) < 3) for clinical in args.column[1:]]):
            input_data = input_data.loc[~(input_data["Stage"] == stage)]
    print(input_data)

    histology &= set(map(step00.get_patient, list(input_data.index)))
    clinical_data = clinical_data.loc[sorted(histology), :]
    print(clinical_data)

    tmp = set(input_data["Stage"])
    order = list(filter(lambda x: x in tmp, step00.long_sample_type_list))
    palette = list(map(lambda x: step00.stage_color_code[x], order))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    for cell in tqdm.tqdm(cells):
        title, tool = cell.split("_")

        try:
            stat, p = scipy.stats.kruskal(*[input_data.loc[(input_data["Stage"] == stage) & (input_data[args.column[0]] == clinical), cell] for clinical in args.column[1:] for stage in order])
        except ValueError:
            _, p = 0.0, 1.0

        if p > 0.05:
            continue

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        seaborn.violinplot(data=input_data, x=args.column[0], order=args.column[1:], y=cell, hue="Stage", hue_order=order, palette=palette, cut=1, linewidth=5, ax=ax)
        statannotations.Annotator.Annotator(ax, [((clinical, a), (clinical, b)) for a, b in itertools.combinations(order, r=2) for clinical in args.column[1:]] + [((a, stage), (b, stage)) for a, b in zip(args.column[1:], args.column[2:]) for stage in order], data=input_data, x=args.column[0], order=args.column[1:], y=cell, hue="Stage", hue_order=order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

        matplotlib.pyplot.title(f"Kruskal-Wallis p={p:.3f}")
        matplotlib.pyplot.ylabel(f"{title} from {tool}")
        matplotlib.pyplot.tight_layout()

        figures.append(cell.replace(" ", "").replace("(", "").replace(")", "").replace("/", "_") + ".pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(figures):
            tar.add(f, arcname=f)
