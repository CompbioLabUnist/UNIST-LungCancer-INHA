"""
draw_clustermap_deconvolution_clinical.py: draw clustermap plot with deconvolution results with clinical data
"""
import argparse
import tarfile
import matplotlib
import matplotlib.colors
import matplotlib.patches
import matplotlib.pyplot
import pandas
import seaborn
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Deconvolution result TSV file (not necessarily TSV)", type=str)
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--column", help="Clinical data column name", type=str, nargs="+", default=["Recurrence", "NO", "YES"])

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0).T
    input_data = input_data.loc[:, list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.columns)))]
    input_data = input_data.loc[list(filter(lambda x: sum(input_data.loc[x, :]) > 0.0, list(input_data.index))), :]
    print(input_data)

    clinical_color_dict = dict(zip(args.column[1:], matplotlib.colors.TABLEAU_COLORS))
    print(clinical_color_dict)

    color_data = pandas.DataFrame(index=input_data.columns)
    color_data["Stage"] = list(map(step00.get_color_by_type, list(color_data.index)))
    color_data[args.column[0]] = list(map(lambda x: clinical_color_dict[clinical_data.loc[step00.get_patient(x), args.column[0]]], list(color_data.index)))
    print(color_data)

    figures = list()

    g = seaborn.clustermap(data=input_data, figsize=(32, 18), row_cluster=True, col_cluster=True, cbar_pos=(-0.04, 0.2, 0.02, 0.6), col_colors=color_data, xticklabels=False, yticklabels=True, square=False, z_score=0, cmap="coolwarm", center=0, robust=True)

    g.ax_heatmap.set_xlabel(f"{input_data.shape[1]} samples")
    g.ax_heatmap.set_ylabel(f"{input_data.shape[0]} cell types")

    matplotlib.pyplot.legend([matplotlib.patches.Patch(facecolor=clinical_color_dict[x]) for x in args.column[1:]], args.column[1:], title=args.column[0], bbox_to_anchor=(0, 1), bbox_transform=matplotlib.pyplot.gcf().transFigure)

    figures.append("All.pdf")
    g.savefig(figures[-1])
    matplotlib.pyplot.close(g.fig)

    stage_set = set(map(step00.get_long_sample_type, list(input_data.columns)))
    stage_list = list(filter(lambda x: x in stage_set, step00.long_sample_type_list))

    for stage in tqdm.tqdm(stage_list):
        selected_samples = list(filter(lambda x: step00.get_long_sample_type(x) == stage, list(input_data.columns)))

        selected_data = input_data.loc[:, selected_samples]
        selected_data = selected_data.loc[list(filter(lambda x: sum(selected_data.loc[x, :]) > 0.0, list(selected_data.index))), :]

        if selected_data.shape[1] < 2:
            continue

        g = seaborn.clustermap(data=selected_data, figsize=(32, 18), row_cluster=True, col_cluster=True, cbar_pos=(-0.04, 0.2, 0.02, 0.6), col_colors=color_data.loc[selected_samples, args.column[0]], xticklabels=False, yticklabels=True, square=False, z_score=0, cmap="coolwarm", center=0, robust=True)

        g.ax_heatmap.set_xlabel(f"{selected_data.shape[1]} {stage} samples")
        g.ax_heatmap.set_ylabel(f"{selected_data.shape[0]} cell types")

        matplotlib.pyplot.legend([matplotlib.patches.Patch(facecolor=clinical_color_dict[x]) for x in args.column[1:]], args.column[1:], title=args.column[0], bbox_to_anchor=(0, 1), bbox_transform=matplotlib.pyplot.gcf().transFigure)

        figures.append(f"{stage}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
