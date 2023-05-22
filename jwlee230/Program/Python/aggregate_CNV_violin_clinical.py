"""
aggregate_CNV_violin_clinical.py: Violin plot of CNV data for PRE-PRI comparing over chromosomes by clinical data
"""
import argparse
import collections
import itertools
import multiprocessing
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

band_data = pandas.DataFrame()


def query_band(chromosome: str, start: int, end: int) -> typing.List[str]:
    answer = list()

    answer += list(band_data.loc[(band_data["chrom"] == chromosome) & (band_data["chrom_start"] <= start) & (start <= band_data["chrom_end"]), "arm"])
    answer += list(band_data.loc[(band_data["chrom"] == chromosome) & (start <= band_data["chrom_start"]) & (band_data["chrom_end"] <= end), "arm"])
    answer += list(band_data.loc[(band_data["chrom"] == chromosome) & (band_data["chrom_start"] <= end) & (end <= band_data["chrom_end"]), "arm"])

    return sorted(set(answer))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="CNV segment.tsv file", type=str)
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("band", help="Chromosome band txt file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--watching", help="Watching column name", type=str, required=True)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)
    parser.add_argument("--compare", help="Comparison grouping (type, control, case)", type=str, nargs=3, default=["Recurrence", "NO", "YES"])

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
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
        control_patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC") & (clinical_data[args.compare[0]] == args.compare[1])].index)
        case_patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC") & (clinical_data[args.compare[0]] == args.compare[2])].index)
    elif args.ADC:
        control_patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC") & (clinical_data[args.compare[0]] == args.compare[1])].index)
        case_patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC") & (clinical_data[args.compare[0]] == args.compare[2])].index)
    else:
        raise Exception("Something went wrong!!")
    patients = set(list(control_patients) + list(case_patients))
    print(sorted(control_patients))
    print(sorted(case_patients))

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data = input_data.loc[(input_data["Patient"].isin(patients))]
    print(input_data)

    sample_list = sorted(set(input_data["Sample"]), key=step00.sorting_by_type)

    with multiprocessing.Pool(args.cpus) as pool:
        input_data["arm"] = pool.starmap(query_band, input_data[["chromosome", "start", "end"]].itertuples(index=False, name=None))
    input_data["length"] = input_data["end"] - input_data["start"] + 1
    input_data = input_data.explode(column="arm", ignore_index=True)
    input_data["chrom"] = list(map(lambda x: "-".join(x), input_data[["chromosome", "arm"]].itertuples(index=False, name=None)))
    print(input_data)

    chromosome_list = list(map(lambda x: "-".join(x), itertools.product(list(filter(lambda x: x in set(input_data["chromosome"]), step00.chromosome_full_list)), ["p", "q"])))
    print(chromosome_list)

    output_data = pandas.DataFrame(data=itertools.product(sample_list, chromosome_list, [1]), columns=["Sample", "Chromosome", args.watching])
    for sample, chromosome in tqdm.tqdm(list(itertools.product(sample_list, chromosome_list))):
        tmp_data = input_data.loc[(input_data["Sample"] == sample) & (input_data["chrom"] == chromosome)]

        if tmp_data.empty:
            continue

        output_data.loc[(output_data["Sample"] == sample) & (output_data["Chromosome"] == chromosome), args.watching] = numpy.average(tmp_data[args.watching], weights=tmp_data["length"])

    output_data["Stage"] = list(map(step00.get_long_sample_type, output_data["Sample"]))
    output_data[args.compare[0]] = list(map(lambda x: args.compare[1] if (step00.get_patient(x) in control_patients) else args.compare[2], output_data["Sample"]))
    print(output_data)

    stage_set = collections.Counter(list(map(step00.get_long_sample_type, sample_list)))
    stage_list = list(filter(lambda x: stage_set[x] > 3, step00.long_sample_type_list))
    print(stage_set)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    ncols = 8
    additional_row = (1 if len(chromosome_list) % ncols else 0)

    fig, axs = matplotlib.pyplot.subplots(ncols=ncols, nrows=len(chromosome_list) // ncols + additional_row, figsize=(ncols * 4 * len(stage_list), 12 * (len(chromosome_list) // ncols + additional_row)))

    for i, chromosome in tqdm.tqdm(list(enumerate(chromosome_list))):
        drawing_data = output_data.loc[(output_data["Chromosome"] == chromosome)]
        drawing_stage_list = list(filter(lambda x: (x in set(drawing_data.loc[(drawing_data[args.compare[0]] == args.compare[1]), "Stage"])) and (x in set(drawing_data.loc[(drawing_data[args.compare[0]] == args.compare[2]), "Stage"])), stage_list))

        pairs = list()
        for stage in drawing_stage_list:
            p = scipy.stats.mannwhitneyu(drawing_data.loc[(drawing_data["Stage"] == stage) & (drawing_data[args.compare[0]] == args.compare[1]), args.watching], drawing_data.loc[(drawing_data["Stage"] == stage) & (drawing_data[args.compare[0]] == args.compare[2]), args.watching])[1]
            if p < 0.05:
                pairs.append(((stage, args.compare[1]), (stage, args.compare[2])))

        try:
            stat, p = scipy.stats.kruskal(*[drawing_data.loc[(drawing_data["Stage"] == stage) & (drawing_data[args.compare[0]] == clinical), args.watching] for stage, clinical in itertools.product(drawing_stage_list, args.compare[1:])])
        except ValueError:
            stat, p = 0.0, 1.0

        seaborn.violinplot(data=drawing_data, x="Stage", y=args.watching, order=drawing_stage_list, hue=args.compare[0], hue_order=args.compare[1:], inner="box", cut=1, ax=axs[i // ncols][i % ncols])
        if pairs:
            statannotations.Annotator.Annotator(axs[i // ncols][i % ncols], pairs, data=drawing_data, x="Stage", y=args.watching, order=drawing_stage_list, hue=args.compare[0], hue_order=args.compare[1:]).configure(test="Mann-Whitney", text_format="star", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()

        axs[i // ncols][i % ncols].set_title(f"{chromosome}: K.W. p={p:.3f}")
        axs[i // ncols][i % ncols].set_xlabel("")
        axs[i // ncols][i % ncols].set_ylabel("Ratio")
        axs[i // ncols][i % ncols].legend(loc="lower left")

    matplotlib.pyplot.tight_layout()
    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
