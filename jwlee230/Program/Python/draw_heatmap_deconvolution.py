"""
draw_heatmap_deconvolution.py: draw heatmap plot with deconvolution results
"""
import argparse
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Deconvolution result TSV file (not necessarily TSV)", type=str)
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

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
    input_data = input_data.loc[sorted(input_data.index), sorted(filter(lambda x: step00.get_patient(x) in patients, list(input_data.columns)), key=step00.sorting_by_type)]
    print(input_data)

    fig, ax = matplotlib.pyplot.subplots(figsize=(40, 18))

    seaborn.heatmap(data=input_data, xticklabels=True, yticklabels=True, square=False, ax=ax, cmap="Reds", vmin=0, vmax=1)

    matplotlib.pyplot.xticks(fontsize="xx-small")
    matplotlib.pyplot.yticks(fontsize="xx-small")
    matplotlib.pyplot.xlabel(f"{input_data.shape[1]} samples")
    matplotlib.pyplot.ylabel(f"{input_data.shape[0]} cell types")

    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
