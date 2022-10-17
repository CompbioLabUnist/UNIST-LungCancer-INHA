"""
draw_clustermap_deconvolution_MSP.py: draw clustermap plot from deconvolution result with Mutation Shared Proportion
"""
import argparse
import itertools
import tarfile
import matplotlib
import matplotlib.colors
import matplotlib.patches
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Deconvolution result TSV file (not necessarily TSV)", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
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
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    print(clinical_data)

    patients = set(clinical_data.index)
    print(patients)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0).T
    input_data = input_data.loc[:, list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.columns)))]
    input_data = input_data.loc[list(filter(lambda x: sum(input_data.loc[x, :]) > 0.0, list(input_data.index))), :]
    print(input_data)

    color_data = pandas.DataFrame(index=input_data.columns)
    color_data["Stage"] = list(map(step00.get_color_by_type, list(color_data.index)))

    for MSP in tqdm.tqdm(step00.sharing_columns):
        if args.median:
            threshold = numpy.median(clinical_data[MSP])
        elif args.mean:
            threshold = numpy.mean(clinical_data[MSP])
        else:
            raise Exception("Something went wrong!!")

        color_data[MSP] = list(map(lambda x: "tab:blue" if (clinical_data.loc[step00.get_patient(x), MSP] < threshold) else "tab:red", list(color_data.index)))

    print(color_data)

    figures = list()
    for MSP in tqdm.tqdm(step00.sharing_columns):
        g = seaborn.clustermap(data=input_data, figsize=(input_data.shape[1], input_data.shape[0]), row_cluster=True, col_cluster=True, cbar_pos=(-0.04, 0.2, 0.02, 0.6), col_colors=color_data[["Stage", MSP]], xticklabels=False, yticklabels=True, square=False, z_score=0, cmap="coolwarm", center=0, robust=True)

        g.ax_heatmap.set_xlabel(f"{input_data.shape[1]} samples")
        g.ax_heatmap.set_ylabel(f"{input_data.shape[0]} cell types")

        matplotlib.pyplot.legend([matplotlib.patches.Patch(facecolor=x) for x in ["tab:blue", "tab:red"]], ["Lower", "Higher"], title=MSP, bbox_to_anchor=(0, 1), bbox_transform=matplotlib.pyplot.gcf().transFigure)

        figures.append(f"{MSP}-All.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

    stage_set = set(map(step00.get_long_sample_type, list(input_data.columns)))
    stage_list = list(filter(lambda x: x in stage_set, step00.long_sample_type_list))

    for stage, MSP in tqdm.tqdm(list(itertools.product(stage_list, step00.sharing_columns))):
        selected_samples = list(filter(lambda x: step00.get_long_sample_type(x) == stage, list(input_data.columns)))

        selected_data = input_data.loc[:, selected_samples]
        selected_data = selected_data.loc[list(filter(lambda x: sum(selected_data.loc[x, :]) > 0.0, list(selected_data.index))), :]

        if selected_data.shape[1] < 2:
            continue

        g = seaborn.clustermap(data=selected_data, figsize=(selected_data.shape[1], selected_data.shape[0]), row_cluster=True, col_cluster=True, cbar_pos=(-0.04, 0.2, 0.02, 0.6), col_colors=color_data.loc[selected_samples, MSP], xticklabels=False, yticklabels=True, square=False, z_score=0, cmap="coolwarm", center=0, robust=True)

        g.ax_heatmap.set_xlabel(f"{selected_data.shape[1]} {stage} samples")
        g.ax_heatmap.set_ylabel(f"{selected_data.shape[0]} cell types")

        matplotlib.pyplot.legend([matplotlib.patches.Patch(facecolor=x) for x in ["tab:blue", "tab:red"]], ["Lower", "Higher"], title=MSP, bbox_to_anchor=(0, 1), bbox_transform=matplotlib.pyplot.gcf().transFigure)

        figures.append(f"{MSP}-{stage}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
