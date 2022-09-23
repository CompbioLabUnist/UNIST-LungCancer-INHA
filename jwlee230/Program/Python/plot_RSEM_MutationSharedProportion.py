"""
plot_RSEM_MutationSharedProportion.py: Plot RSEM gene expression with Mutation Shared Proportion
"""
import argparse
import itertools
import multiprocessing
import tarfile
import typing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import statannotations.Annotator
import tqdm
import step00

input_data = pandas.DataFrame()
gene_list = list()


def pearson(sample: str) -> typing.Tuple[str, float, float]:
    drawing_data = input_data.loc[[sample, step00.get_paired_primary(sample)], gene_list].T
    fig_name = f"Joint_{sample}_{step00.get_paired_primary(sample)}.pdf"

    r, p = scipy.stats.pearsonr(drawing_data[sample], drawing_data[step00.get_paired_primary(sample)])

    g = seaborn.jointplot(data=drawing_data, x=sample, y=step00.get_paired_primary(sample), kind="reg", height=24, dropna=True)
    g.fig.text(0.5, 0.75, "r={0:.3f}, p={1:.3f}".format(r, p), color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")
    g.set_axis_labels(f"{sample} ({step00.get_long_sample_type(sample)})", f"{step00.get_paired_primary(sample)} (Primary)")

    g.savefig(fig_name)
    matplotlib.pyplot.close(g.fig)

    return fig_name, r, p


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="RSEM expression TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

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

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0).T
    gene_list = list(input_data.columns)
    print(input_data)

    with multiprocessing.Pool(args.cpus) as pool:
        input_data["Patient"] = pool.map(step00.get_patient, list(input_data.index))
        input_data["Stage"] = pool.map(step00.get_long_sample_type, list(input_data.index))
    input_data = input_data.loc[(input_data["Patient"].isin(patients))]
    print(input_data)

    patients &= set(input_data["Patient"])
    print(patients)

    sample_list = list(input_data.index)
    precancer_list = list(filter(lambda x: (step00.get_long_sample_type(x) == "CIS+AIS") and (step00.get_paired_primary(x) in sample_list), sample_list))
    precancer_patient_list = list(map(step00.get_patient, precancer_list))
    print(precancer_list)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    output_data = pandas.DataFrame(index=precancer_patient_list)
    with multiprocessing.Pool(args.cpus) as pool:
        pool.map(pearson, precancer_list)
        output_data[["Joint", "r", "p"]] = pool.map(pearson, precancer_list)
    print(output_data)

    figures = list(output_data["Joint"])
    del output_data["Joint"]

    with multiprocessing.Pool(args.cpus) as pool:
        output_data["r"] = pool.map(float, output_data["r"])
        output_data["p"] = pool.map(float, output_data["p"])

    output_data["-log10(p)"] = -1 * numpy.log10(output_data["p"].to_numpy())
    print(output_data)

    for MSP in tqdm.tqdm(step00.sharing_columns):
        if args.median:
            threshold = numpy.median(clinical_data[MSP])
        elif args.mean:
            threshold = numpy.mean(clinical_data[MSP])
        else:
            raise ValueError("Something went wrong!!")

        output_data[MSP] = list(map(lambda x: "Lower" if (clinical_data.loc[x, MSP] < threshold) else "Higher", precancer_patient_list))
    print(output_data)

    for MSP in tqdm.tqdm(step00.sharing_columns):
        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.scatterplot(data=output_data, x="r", y="-log10(p)", hue=MSP, hue_order=["Lower", "Higher"], palette={"Lower": "tab:blue", "Higher": "tab:red"}, s=1000, legend="full", ax=ax)

        matplotlib.pyplot.tight_layout()

        figures.append(f"Scatter_{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
