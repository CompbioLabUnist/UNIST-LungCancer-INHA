"""
draw_Tcell_density_violin_deconvolution_MSP.py: draw T-cell density violin plot for deconvolution result with Mutation Shared Proportion
"""
import argparse
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

    parser.add_argument("input", help="Deconvolution result TSV file (not necessarily TSV)", type=str)
    parser.add_argument("clinical", help="Clinical data w/ mutation shared proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_threshold = parser.add_mutually_exclusive_group(required=True)
    group_threshold.add_argument("--median", help="Use median threshold", action="store_true", default=False)
    group_threshold.add_argument("--mean", help="Use mean threshold", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clincal_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clincal_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")

    patients = set(clinical_data.index)
    print(patients)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)

    input_data = input_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index))), step00.Tcell_list]
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
    MSP_order = ["Lower", "Higher"]

    figures = list()

    for MSP in tqdm.tqdm(step00.sharing_columns):
        if args.median:
            threshold = numpy.median(clinical_data[MSP])
        elif args.mean:
            threshold = numpy.mean(clinical_data[MSP])
        else:
            raise Exception("Something went wrong!!")

        input_data[MSP] = list(map(lambda x: "Lower" if (clinical_data.loc[step00.get_patient(x), MSP] < threshold) else "Higher", list(input_data.index)))

        try:
            stat, p = scipy.stats.kruskal(*[input_data.loc[(input_data["Stage"] == stage) & (input_data[MSP] == clinical), "T cell density"] for clinical in MSP_order for stage in order])
        except ValueError:
            stat, p = 0.0, 1.0

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        seaborn.violinplot(data=input_data, x=MSP, y="T cell density", order=MSP_order, hue="Stage", hue_order=order, palette=palette, cut=1, linewidth=5, ax=ax)
        statannotations.Annotator.Annotator(ax, [((clinical, a), (clinical, b)) for a, b in zip(order, order[1:]) for clinical in MSP_order] + [((a, stage), (b, stage)) for a, b in zip(MSP_order, MSP_order[1:]) for stage in order], data=input_data, x=MSP, y="T cell density", order=MSP_order, hue="Stage", hue_order=order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

        matplotlib.pyplot.title(f"Kruskal-Wallis p={p:.3f}")
        matplotlib.pyplot.tight_layout()

        figures.append(f"{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
