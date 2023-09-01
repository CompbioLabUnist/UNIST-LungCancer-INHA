"""
draw_bar_plot_deconvolution_MSP.py: draw bar plot from deconvolution results with MSP
"""
import argparse
import itertools
import tarfile
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import numpy
import pandas
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Deconvolution result TSV file (not necessarily TSV)", type=str)
    parser.add_argument("clinical", help="Clinical data w/ Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    else:
        raise Exception("Something went wrong!!")
    print(clinical_data)

    patients = set(clinical_data.index)
    print(patients)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data = input_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index))), :]
    print(input_data)

    cells = sorted(input_data.columns)
    cell_palette = dict(zip(cells, itertools.cycle(matplotlib.colors.XKCD_COLORS)))

    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    for MSP in tqdm.tqdm(step00.sharing_columns):
        input_data[MSP] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), MSP], list(input_data.index)))
    print(input_data)

    figures = list()
    for stage, MSP in tqdm.tqdm(list(itertools.product(step00.long_sample_type_list, step00.sharing_columns))):
        drawing_data = input_data.loc[(input_data["Stage"] == stage)].sort_values(MSP)

        if drawing_data.empty:
            continue

        cells = sorted(cells, key=lambda x: numpy.mean(drawing_data[x]), reverse=True)
        drawing_data = drawing_data.loc[:, cells]

        fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

        for i, index in enumerate(list(drawing_data.index)):
            for j, cell in enumerate(cells):
                labeling = (i == 0) and (j < 5)
                matplotlib.pyplot.bar(x=i, height=drawing_data.iloc[i, j], bottom=sum(drawing_data.iloc[i, :j]), color=cell_palette[cell], label=cell if labeling else None, edgecolor=None, linewidth=0)

        matplotlib.pyplot.xticks([])
        matplotlib.pyplot.xlabel(f"{len(drawing_data)} {stage} samples ordered by {MSP}")
        matplotlib.pyplot.ylabel("Cell type proportions")
        matplotlib.pyplot.ylim(0, 1)
        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.legend(loc="upper right", fontsize="xx-small")
        matplotlib.pyplot.tight_layout()

        figures.append(f"{stage}-{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
