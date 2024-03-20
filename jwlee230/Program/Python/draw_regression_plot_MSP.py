"""
draw_regression_plot_MSP.py: draw regression plot for MSP
"""
import argparse
import itertools
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import scipy.stats
import seaborn
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input expression TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(len(patients), patients)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0).T
    input_data = input_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index))), :]
    gene_list = list(input_data.columns)
    sample_list = list(input_data.index)
    print(input_data)

    for MSP in tqdm.tqdm(step00.sharing_columns):
        input_data[MSP] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), MSP], list(input_data.index)))
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    for MSP, gene in tqdm.tqdm(list(itertools.product(step00.sharing_columns[:1], ["H2AC4"]))):
        precancer_list = list(filter(lambda x: x in sample_list, list(clinical_data[f"{MSP}-sample"])))
        tmp_data = input_data.loc[precancer_list, [gene, MSP]]

        print(gene, scipy.stats.linregress(tmp_data[MSP], tmp_data[gene]))

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.regplot(data=tmp_data, x=MSP, y=gene, scatter=True, fit_reg=True, color="tab:pink", ax=ax)

        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.tight_layout()

        figures.append(f"{gene}-{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
