"""
draw_clustermap_RSEM.py: draw clustermap upon RSEM DEG data
"""
import argparse
import matplotlib
import matplotlib.pyplot
import pandas
import scipy
import seaborn
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="RSEM result TSV file", type=str)
    parser.add_argument("clinical", help="Clinical CSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col="gene_name")
    print(input_data)

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    patients = list(input_data.columns)
    SQCs = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    ADCs = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    patients = list(filter(lambda x: (step00.get_patient(x) in SQCs) or (step00.get_patient(x) in ADCs), patients))
    print(patients)

    input_data = input_data.loc[:, patients]
    for index in list(input_data.index):
        input_data.loc[index, :] = scipy.stats.zscore(input_data.loc[index, :])
    input_data.dropna(axis="index", how="any", inplace=True)

    input_data["std"] = list(input_data.std(axis="columns"))
    input_data.sort_values(by="std", ascending=False, inplace=True)
    input_data = input_data.iloc[:len(input_data) // 100, :]
    del input_data["std"]
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    g = seaborn.clustermap(data=input_data, figsize=(input_data.shape[1], input_data.shape[0]), row_cluster=True, col_cluster=True, cbar_pos=None, col_colors=list(map(lambda x: "tab:red" if step00.get_patient(x) in SQCs else "tab:blue", input_data.columns)), tree_kws={"linewidths": 2.0}, xticklabels=True, yticklabels=True, square=False, cmap="bwr")
    g.ax_heatmap.set_ylabel("")

    g.savefig(args.output)
