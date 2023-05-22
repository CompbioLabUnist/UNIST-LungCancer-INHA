"""
plot_CNV_MutationSharedProportion.py: violin plot & joint plot Copy Number Variation with Mutation Shared Proportion
"""
import argparse
import itertools
import multiprocessing
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

    parser.add_argument("input", help="CNV ploidy.tsv file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)
    parser.add_argument("--threshold", help="Threshold for gain/loss", type=float, default=0.2)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_strategy = parser.add_mutually_exclusive_group(required=True)
    group_strategy.add_argument("--median", help="Median division", action="store_true", default=False)
    group_strategy.add_argument("--mean", help="Mean division", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.threshold < 1):
        raise ValueError("Threshold must be (0, 1)")

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data["Patient"] = list(map(step00.get_patient, list(input_data.index)))
    input_data = input_data.loc[(input_data["Patient"].isin(patients))]
    print(input_data)

    sample_list = sorted(set(input_data.index), key=step00.sorting_by_type)
    print(sample_list)

    output_data = pandas.DataFrame(data=sample_list, columns=["Sample"])
    with multiprocessing.Pool(args.cpus) as pool:
        output_data["Patient"] = pool.map(step00.get_patient, output_data["Sample"])
        output_data["Stage"] = pool.map(step00.get_long_sample_type, output_data["Sample"])
        output_data["Ploidy"] = list(map(lambda x: input_data.loc[x, "Ploidy"], output_data["Sample"]))

    sample_list = sorted(set(output_data["Sample"]), key=step00.sorting_by_type)
    print(sample_list)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    MSP_order = ["Lower", "Higher"]

    for MSP in tqdm.tqdm(step00.sharing_columns):
        if args.median:
            cutting = numpy.median(clinical_data[MSP])
        elif args.mean:
            cutting = numpy.mean(clinical_data[MSP])
        else:
            raise Exception("Something went wrong!!")

        output_data[MSP] = list(map(lambda x: "Lower" if (clinical_data.loc[x, MSP] < cutting) else "Higher", output_data["Patient"]))

        stage_list = list(filter(lambda x: all([(x in set(output_data.loc[(output_data[MSP] == x1), "Stage"])) for x1 in MSP_order]), step00.long_sample_type_list))
        palette = dict([(stage, step00.stage_color_code[stage]) for stage in stage_list])

        compare_list = list()
        for (x1, s1), (x2, s2) in [(("Lower", stage), ("Higher", stage)) for stage in stage_list] + [((x, a), (x, b)) for x in MSP_order for a, b in itertools.combinations(stage_list, r=2)]:
            stat, p = scipy.stats.mannwhitneyu(output_data.loc[(output_data[MSP] == x1) & (output_data["Stage"] == s1), "Ploidy"], output_data.loc[(output_data[MSP] == x2) & (output_data["Stage"] == s2), "Ploidy"])
            if p < 0.05:
                compare_list.append(((x1, s1), (x2, s2)))

        try:
            stat, p = scipy.stats.kruskal(*[output_data.loc[(output_data[MSP] == x) & (output_data["Stage"] == s), "Ploidy"] for x, s in itertools.product(MSP_order, stage_list)])
        except ValueError:
            p = 1.0

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.violinplot(data=output_data, x=MSP, order=MSP_order, y="Ploidy", hue="Stage", hue_order=stage_list, palette=palette, inner="box", cut=1, ax=ax)
        if compare_list:
            statannotations.Annotator.Annotator(ax, compare_list, data=output_data, x=MSP, order=["Lower", "Higher"], y="Ploidy", hue="Stage", hue_order=stage_list).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()

        matplotlib.pyplot.title(f"Kruskal-Wallis p={p:.3f}")
        matplotlib.pyplot.tight_layout()

        figures.append(f"Violin_{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

        del output_data[MSP]

    stage_list = list(filter(lambda x: output_data.loc[(output_data["Stage"] == x)].shape[0] > 3, step00.long_sample_type_list))
    palette = dict([(stage, step00.stage_color_code[stage]) for stage in stage_list])

    output_data = output_data.loc[(output_data["Stage"].isin(stage_list))]

    for MSP in tqdm.tqdm(step00.sharing_columns):
        output_data[MSP] = list(map(lambda x: clinical_data.loc[x, MSP], output_data["Patient"]))

        g = seaborn.jointplot(data=output_data, x=MSP, y="Ploidy", hue="Stage", hue_order=stage_list, palette=palette, height=24, ratio=5, kind="scatter")
        figures.append(f"Joint_All_{MSP}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

    for stage, MSP in tqdm.tqdm(list(itertools.product(stage_list, step00.sharing_columns))):
        tmp_data = output_data.loc[(output_data["Stage"] == stage)]

        r, p = scipy.stats.pearsonr(tmp_data[MSP], tmp_data["Ploidy"])

        g = seaborn.jointplot(data=tmp_data, x=MSP, y="Ploidy", color=palette[stage], height=24, ratio=5, kind="reg")
        g.fig.text(0.5, 0.75, "r={0:.3f}, p={1:.3f}".format(r, p), color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")
        figures.append(f"Joint_{stage}_{MSP}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
