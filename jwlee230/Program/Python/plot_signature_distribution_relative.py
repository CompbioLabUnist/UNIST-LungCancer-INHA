"""
plot_signature_distribution_relative.py: plot signature distribution in relative
"""
import argparse
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Signature TSV file (not necessarily TSV)", type=str)
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".csv"):
        raise ValueError("Clinical data must end with .csv!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col="Samples")
    signatures = list(input_data.columns)
    input_data["Total"] = input_data.sum(axis="columns")
    for index in list(input_data.index):
        input_data.loc[index, :] = input_data.loc[index, :] / input_data.loc[index, "Total"]
    print(input_data)

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(sorted(patients))

    input_data = input_data.loc[sorted(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index)), key=step00.sorting_by_type), :]
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    order = list(filter(lambda x: x in set(input_data["Stage"]), step00.long_sample_type_list))
    print(input_data)
    print(order)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    g = seaborn.pairplot(data=input_data, hue="Stage", hue_order=order, palette=step00.stage_color_code, x_vars=signatures, y_vars=reversed(signatures), kind="scatter", diag_kind="kde", height=8, aspect=1)

    g.tight_layout()

    g.savefig(args.output)
    matplotlib.pyplot.close(g.fig)
