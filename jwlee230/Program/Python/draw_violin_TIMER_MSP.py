"""
draw_violin_TIMER_MSP.py: Draw violin plot from TIMER withi Mutation Shared Proportion
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="TIMER result CSV file", type=str)
    parser.add_argument("clinical", help="Clinical data data w/ Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

    group_threshold = parser.add_mutually_exclusive_group(required=True)
    group_threshold.add_argument("--median", help="Use median threshold", action="store_true", default=False)
    group_threshold.add_argument("--mean", help="Use mean threshold", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".csv"):
        raise ValueError("Input must end with .CSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    histology = set(clinical_data.index)
    print(histology)

    input_data = pandas.read_csv(args.input, index_col=0).T
    cells = list(input_data.columns)
    print(input_data)

    histology &= set(map(step00.get_patient, list(input_data.index)))
    clinical_data = clinical_data.loc[sorted(histology), :]
    print(clinical_data)

    input_data = input_data.loc[list(filter(lambda x: step00.get_patient(x) in histology, list(input_data.index))), :]
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    print(input_data)

    for MSP in tqdm.tqdm(step00.sharing_columns):
        if args.median:
            threshold = numpy.median(clinical_data[MSP])
        elif args.mean:
            threshold = numpy.mean(clinical_data[MSP])
        else:
            raise Exception("Something went wrong!!")

        input_data[MSP] = list(map(lambda x: "Lower" if (clinical_data.loc[step00.get_patient(x), MSP] < threshold) else "Higher", list(input_data.index)))
    print(input_data)

    MSP_order = ["Lower", "Higher"]

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    for cell, MSP in tqdm.tqdm(list(itertools.product(cells, step00.sharing_columns))):
        title, tool = cell.split("_")
        tmp = set.intersection(*[set(input_data.loc[(input_data[MSP] == clinical), "Stage"]) for clinical in MSP_order])
        order = list(filter(lambda x: x in tmp, step00.long_sample_type_list))
        palette = list(map(lambda x: step00.stage_color_code[x], order))

        try:
            stat, p = scipy.stats.kruskal(*[input_data.loc[(input_data["Stage"] == stage) & (input_data[MSP] == clinical), cell] for clinical in MSP_order for stage in order])
        except ValueError:
            _, p = 0.0, 1.0

        if p > 0.05:
            continue

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        seaborn.violinplot(data=input_data, x=MSP, order=MSP_order, y=cell, hue="Stage", hue_order=order, palette=palette, cut=1, linewidth=5, ax=ax)
        statannotations.Annotator.Annotator(ax, [((clinical, a), (clinical, b)) for a, b in itertools.combinations(order, r=2) for clinical in MSP_order] + [((a, stage), (b, stage)) for a, b in zip(MSP_order, MSP_order[1:]) for stage in order], data=input_data, x=MSP, order=MSP_order, y=cell, hue="Stage", hue_order=order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

        matplotlib.pyplot.title(f"Kruskal-Wallis p={p:.3f}")
        matplotlib.pyplot.ylabel(f"{title} from {tool}")
        matplotlib.pyplot.tight_layout()

        figures.append(f"{MSP}-{tool}-{title}.pdf".replace("/", "_"))
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(figures):
            tar.add(f, arcname=f)
