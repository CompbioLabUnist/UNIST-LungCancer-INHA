"""
aggregate_sequenza_violin.py: Violin plot Sequenza data for PRE-PRI comparing over chromosomes
"""
import argparse
import itertools
import multiprocessing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import statannotations.Annotator
import tqdm
import step00

watching = "depth.ratio"


def get_data(file_name: str) -> pandas.DataFrame:
    data = pandas.read_csv(file_name, sep="\t", usecols=["chromosome", "start.pos", "end.pos", watching]).dropna(axis="index")
    data["sample"] = file_name.split("/")[-2]
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Sequenza output segments.txt file(s)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".txt"), args.input)):
        raise ValueError("INPUT must end with .TXT!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
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
    print(patients)

    args.input = list(filter(lambda x: step00.get_patient(x.split("/")[-2]) in patients, args.input))
    sample_list = list(map(lambda x: x.split("/")[-2], args.input))
    print(sample_list)

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="index", copy=False, ignore_index=True, verify_integrity=True)
    input_data["length"] = input_data["end.pos"] - input_data["start.pos"] + 1
    print(input_data)

    chromosome_list = list(filter(lambda x: x in set(input_data["chromosome"]), step00.chromosome_list))
    print(chromosome_list)

    output_data = pandas.DataFrame(data=itertools.product(sample_list, chromosome_list, [0]), columns=["Sample", "Chromosome", watching])
    for sample, chromosome in tqdm.tqdm(itertools.product(sample_list, chromosome_list)):
        tmp_data = input_data.loc[(input_data["sample"] == sample) & (input_data["chromosome"] == chromosome)]
        output_data.loc[(output_data["Sample"] == sample) & (output_data["Chromosome"] == chromosome), watching] = numpy.average(tmp_data[watching], weights=tmp_data["length"])
    output_data["PRE/PRI"] = list(map(step00.get_simple_sample_type, output_data["Sample"]))
    print(output_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, axs = matplotlib.pyplot.subplots(ncols=6, nrows=len(chromosome_list) // 6 + (1 if len(chromosome_list) % 6 else 0), figsize=(64, 48))

    for i, chromosome in tqdm.tqdm(enumerate(chromosome_list)):
        drawing_data = output_data.loc[(output_data["Chromosome"] == chromosome)]

        seaborn.violinplot(data=drawing_data, x="PRE/PRI", y=watching, order=["Precancer", "Primary"], inner="box", ax=axs[i // 6][i % 6])
        statannotations.Annotator.Annotator(axs[i // 6][i % 6], [("Precancer", "Primary")], data=drawing_data, x="PRE/PRI", y=watching, order=["Precancer", "Primary"]).configure(test="Mann-Whitney", text_format="star", loc="inside", verbose=0).apply_and_annotate()

        axs[i // 6][i % 6].set_title(chromosome)

    matplotlib.pyplot.tight_layout()
    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
