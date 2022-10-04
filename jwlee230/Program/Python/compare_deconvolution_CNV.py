"""
compare_deconvolution_CNV.py: compare deconvolution result with CNV result
"""
import argparse
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Deconvolution result TSV file", type=str)
    parser.add_argument("CNV", help="Sequenza output sample.tsv file", type=str)
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Deconvolution file must end with .TSV!!")
    elif not args.CNV.endswith(".tsv"):
        raise ValueError("CNV file must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    CNV_data = pandas.read_csv(args.CNV, sep="\t", index_col=0)
    print(CNV_data)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    deconvolution_cells = sorted(input_data.columns)
    print(input_data)

    output_data = pandas.concat([input_data, CNV_data], axis="columns", join="inner", verify_integrity=True, sort=True)
    print(output_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    stage_list = set(output_data["Stage"])
    order = list(filter(lambda x: x in stage_list, step00.long_sample_type_list))
    palette = list(map(lambda x: step00.stage_color_code[x], order))

    figures = list()
    for cell in tqdm.tqdm(deconvolution_cells):
        if numpy.var(output_data[cell]) == 0:
            continue

        r, p = scipy.stats.pearsonr(output_data["Ploidy"], output_data[cell])

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        g = seaborn.jointplot(data=output_data, x="Ploidy", y=cell, kind="scatter", height=24, ratio=6, hue="Stage", palette=palette, hue_order=order)
        g.fig.text(0.5, 0.75, "r={0:.3f}, p={1:.3f}".format(r, p), color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")
        g.plot_marginals(seaborn.histplot, kde=True, stat="probability", multiple="stack")
        g.set_axis_labels("Ploidy", "{0} proportion".format(cell))

        figures.append(cell.replace(" ", "").replace("(", "").replace(")", "").replace("/", "_") + ".pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close()

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
