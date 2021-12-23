"""
aggregate_PureCN_violin.py: Violin plot PureCN data for PRE-PRI comparing over chromosomes
"""
import argparse
import itertools
import multiprocessing
import typing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import statannotations.Annotator
import tqdm
import step00

watching = "seg.mean"
band_data = pandas.DataFrame()


def get_data(file_name: str) -> pandas.DataFrame:
    data = pandas.read_csv(file_name, sep="\t")
    data["chrom"] = list(map(lambda x: "chrX" if (str(x) == "23") else ("chr" + str(x)), data["chrom"]))
    data["ID"] = step00.get_id(file_name)
    data[watching] = numpy.power(2, data[watching])
    return data


def query_band(chromosome: str, start: int, end: int) -> typing.List[str]:
    answer = list()

    answer += list(band_data.loc[(band_data["chrom"] == chromosome) & (band_data["chrom_start"] <= start) & (start <= band_data["chrom_end"]), "arm"])
    answer += list(band_data.loc[(band_data["chrom"] == chromosome) & (start <= band_data["chrom_start"]) & (band_data["chrom_end"] <= end), "arm"])
    answer += list(band_data.loc[(band_data["chrom"] == chromosome) & (band_data["chrom_start"] <= end) & (end <= band_data["chrom_end"]), "arm"])

    return sorted(set(answer))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PureCN output segments.TSV file(s)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("band", help="Chromosome band txt file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    band_data = step00.get_band_data(args.band)
    print(band_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))
    sample_list = list(map(step00.get_id, args.input))
    print(sample_list)

    stage_set = set(map(step00.get_long_sample_type, args.input))
    stage_list = list(filter(lambda x: x in stage_set, step00.long_sample_type_list))

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="index", copy=False, ignore_index=True, verify_integrity=True)
        input_data["arm"] = pool.starmap(query_band, input_data[["chrom", "loc.start", "loc.end"]].itertuples(index=False, name=None))
    input_data["length"] = input_data["loc.end"] - input_data["loc.start"] + 1
    input_data = input_data.explode(column="arm", ignore_index=True).dropna(axis="index", subset=["arm"])
    input_data["chrom-arm"] = list(map(lambda x: "-".join(x), input_data[["chrom", "arm"]].itertuples(index=False, name=None)))
    print(input_data)

    chromosome_list = list(map(lambda x: "-".join(x), itertools.product(list(filter(lambda x: x in set(input_data["chrom"]), step00.chromosome_full_list)), ["p", "q"])))
    print(chromosome_list)

    output_data = pandas.DataFrame(data=itertools.product(sample_list, chromosome_list, [1]), columns=["Sample", "Chromosome", watching])
    for sample, chromosome in tqdm.tqdm(list(itertools.product(sample_list, chromosome_list))):
        tmp_data = input_data.loc[(input_data["ID"] == sample) & (input_data["chrom-arm"] == chromosome)]

        if tmp_data.empty:
            continue

        output_data.loc[(output_data["Sample"] == sample) & (output_data["Chromosome"] == chromosome), watching] = numpy.average(tmp_data[watching], weights=tmp_data["length"])

    output_data["Stage"] = list(map(step00.get_long_sample_type, output_data["Sample"]))
    print(output_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    ncols = 8
    additional_row = (1 if len(chromosome_list) % ncols else 0)

    fig, axs = matplotlib.pyplot.subplots(ncols=ncols, nrows=len(chromosome_list) // ncols + additional_row, figsize=(ncols * 4 * len(stage_list), 12 * (len(chromosome_list) // ncols + additional_row)))

    for i, chromosome in tqdm.tqdm(list(enumerate(chromosome_list))):
        drawing_data = output_data.loc[(output_data["Chromosome"] == chromosome)]
        drawing_stage_list = list(filter(lambda x: x in set(drawing_data["Stage"]), stage_list))

        seaborn.violinplot(data=drawing_data, x="Stage", y=watching, order=drawing_stage_list, inner="box", ax=axs[i // ncols][i % ncols])
        statannotations.Annotator.Annotator(axs[i // ncols][i % ncols], list(zip(drawing_stage_list, drawing_stage_list[1:])), data=drawing_data, x="Stage", y=watching, order=drawing_stage_list).configure(test="Mann-Whitney", text_format="star", loc="inside", verbose=0).apply_and_annotate()

        axs[i // ncols][i % ncols].set_title(chromosome)
        axs[i // ncols][i % ncols].set_xlabel("")

    matplotlib.pyplot.tight_layout()
    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
