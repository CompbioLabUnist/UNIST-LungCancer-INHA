"""
plot_RSEM_MutationSharedProportion_distance.py: Plot RSEM gene expression with Mutation Shared Proportion with Distance
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
gene_list: typing.List[str] = list()


def L1(sample: str) -> float:
    return sum(list(map(lambda x: abs(x[0] - x[1]), zip(input_data.loc[sample, gene_list], input_data.loc[step00.get_paired_primary(sample), gene_list]))))


def L2(sample: str) -> float:
    return numpy.sqrt(sum(list(map(lambda x: (x[0] - x[1]) ** 2, zip(input_data.loc[sample, gene_list], input_data.loc[step00.get_paired_primary(sample), gene_list])))))


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
        output_data["L1"] = pool.map(L1, precancer_list)
        output_data["L2"] = pool.map(L2, precancer_list)
    print(output_data)

    distances = ["L1", "L2"]
    with multiprocessing.Pool(args.cpus) as pool:
        for distance in tqdm.tqdm(distances):
            output_data[distance] = pool.map(float, output_data[distance])
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

    figures = list()

    for MSP, distance in tqdm.tqdm(list(itertools.product(step00.sharing_columns, distances))):
        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.violinplot(data=output_data, x=MSP, order=["Lower", "Higher"], y=f"{distance}", inner="box", cut=1, ax=ax)
        statannotations.Annotator.Annotator(ax, [("Lower", "Higher")], data=output_data, x=MSP, order=["Lower", "Higher"], y=f"{distance}").configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

        matplotlib.pyplot.tight_layout()

        figures.append(f"Violin_{distance}_{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
