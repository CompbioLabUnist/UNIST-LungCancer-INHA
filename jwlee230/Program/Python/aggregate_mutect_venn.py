"""
aggregate_mutect_venn.py: aggregate mutect MAF files as venn diagram
"""
import argparse
import multiprocessing
import tarfile
import typing
import matplotlib
import matplotlib.pyplot
import seaborn
import pandas
import tqdm
import upsetplot
import step00

patient_samples: typing.Dict[str, typing.Set[str]] = dict()
query_data: typing.Dict[str, typing.Set[str]] = dict()


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


def draw_venn(sample: str) -> str:
    figure_name = f"{sample}.pdf"
    sample_list = sorted(patient_samples[sample], key=step00.sorting_by_type)
    venn_data = dict()

    for sample in sample_list:
        venn_data[sample] = query_data[sample]

    fig = matplotlib.pyplot.figure(figsize=(2 ** len(venn_data) + 40, 24))

    try:
        upsetplot.plot(upsetplot.from_contents(venn_data), fig=fig, show_counts="%d", show_percentages=True, element_size=None)
    except AttributeError:
        pass

    fig.savefig(figure_name, bbox_inches="tight")
    matplotlib.pyplot.close(fig)

    return figure_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(sorted(patients))

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))
    sample_list = list(map(step00.get_id, args.input))

    for sample in tqdm.tqdm(sample_list):
        if step00.get_patient(sample) not in patient_samples:
            patient_samples[step00.get_patient(sample)] = set()
        patient_samples[step00.get_patient(sample)].add(sample)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
        mutect_data["Tumor_Sample_Barcode"] = pool.map(step00.get_id, mutect_data["Tumor_Sample_Barcode"])
    print(mutect_data)

    for sample in tqdm.tqdm(sample_list):
        query_data[sample] = set(mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == sample), step00.sharing_strategy].itertuples(index=False, name=None))

    with multiprocessing.Pool(args.cpus) as pool:
        figures = pool.map(draw_venn, list(patient_samples.keys()))

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
