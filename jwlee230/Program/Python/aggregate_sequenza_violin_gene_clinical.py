"""
aggregate_sequenza_violin_gene_clinical.py: Violin plot Sequenza data for PRE-PRI comparing over genes by clinical data
"""
import argparse
import itertools
import tarfile
import multiprocessing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import statannotations.Annotator
import tqdm
import step00

watching = "depth.ratio"
band_data = pandas.DataFrame()


def get_data(file_name: str) -> pandas.DataFrame:
    data = pandas.read_csv(file_name, sep="\t", usecols=["chromosome", "start.pos", "end.pos", watching]).dropna(axis="index")
    data["sample"] = file_name.split("/")[-2]
    return data


def get_chromosome(location: str) -> str:
    return "chr" + location.split(":")[0]


def get_start(location: str) -> int:
    return int(location.replace("-", ":").split(":")[1])


def get_end(location: str) -> int:
    return int(location.replace("-", ":").split(":")[2])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Sequenza output segments.txt file(s)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("cgc", help="CGC CSV files", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--compare", help="Comparison grouping (type, control, case)", type=str, nargs=3, default=["Recurrence", "NO", "YES"])
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".txt"), args.input)):
        raise ValueError("INPUT must end with .TXT!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.cgc.endswith(".csv"):
        raise ValueError("CGC must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

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

    args.input = list(filter(lambda x: step00.get_patient(x.split("/")[-2]) in patients, args.input))
    sample_list = list(map(lambda x: x.split("/")[-2], args.input))
    print(sample_list)

    stage_set = set(map(step00.get_long_sample_type, sample_list))
    stage_list = list(filter(lambda x: x in stage_set, step00.long_sample_type_list))

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="index", copy=False, ignore_index=True, verify_integrity=True)
    input_data["length"] = input_data["end.pos"] - input_data["start.pos"] + 1
    print(input_data)

    cgc_data = pandas.read_csv(args.cgc)
    cgc_data = cgc_data.loc[~(cgc_data["Genome Location"].str.contains(":-"))]
    with multiprocessing.Pool(args.cpus) as pool:
        cgc_data["Chromosome"] = pool.map(get_chromosome, cgc_data["Genome Location"])
        cgc_data["Start"] = pool.map(get_start, cgc_data["Genome Location"])
        cgc_data["End"] = pool.map(get_end, cgc_data["Genome Location"])
    cgc_data.set_index("Gene Symbol", verify_integrity=True, inplace=True)
    print(cgc_data)

    gene_list = list(cgc_data.index)

    chromosome_list = list(filter(lambda x: x in set(input_data["chromosome"]), step00.chromosome_full_list))
    print(chromosome_list)

    output_data = pandas.DataFrame(data=itertools.product(sample_list, gene_list, [1]), columns=["Sample", "Gene", watching])
    for sample, gene in tqdm.tqdm(list(itertools.product(sample_list, gene_list))):
        chromosome = cgc_data.loc[gene, "Chromosome"]
        start = cgc_data.loc[gene, "Start"]
        end = cgc_data.loc[gene, "End"]
        length = end - start + 1

        a = list()
        weights = list()

        tmp_data = input_data.loc[(input_data["sample"] == sample) & (input_data["chromosome"] == chromosome) & (input_data["start.pos"] <= start) & (end <= input_data["end.pos"]), :]
        for index, row in tmp_data.iterrows():
            a.append(row[watching])
            weights.append(1)

        tmp_data = input_data.loc[(input_data["sample"] == sample) & (input_data["chromosome"] == chromosome) & (start <= input_data["start.pos"]) & (input_data["end.pos"] <= end), :]
        for index, row in tmp_data.iterrows():
            a.append(row[watching])
            weights.append((row["end.pos"] - row["start.pos"] + 1) / length)

        tmp_data = input_data.loc[(input_data["sample"] == sample) & (input_data["chromosome"] == chromosome) & (input_data["start.pos"] <= start) & (start <= input_data["end.pos"]), :]
        for index, row in tmp_data.iterrows():
            a.append(row[watching])
            weights.append((row["end.pos"] - start + 1) / length)

        tmp_data = input_data.loc[(input_data["sample"] == sample) & (input_data["chromosome"] == chromosome) & (input_data["start.pos"] <= end) & (end <= input_data["end.pos"]), :]
        for index, row in tmp_data.iterrows():
            a.append(row[watching])
            weights.append((end - row["start.pos"] + 1) / length)

        if a and weights:
            output_data.loc[(output_data["Sample"] == sample) & (output_data["Gene"] == gene), watching] = numpy.average(a=a, weights=weights)

    output_data["Stage"] = list(map(step00.get_long_sample_type, output_data["Sample"]))
    output_data[args.compare[0]] = list(map(lambda x: args.compare[1] if step00.get_patient(x) in control_patients else args.compare[2], output_data["Sample"]))
    print(output_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()

    for gene in tqdm.tqdm(gene_list):
        fig, ax = matplotlib.pyplot.subplots(figsize=(5 * len(stage_list), 18))

        drawing_data = output_data.loc[(output_data["Gene"] == gene), :]
        drawing_stage_list = list(filter(lambda x: (x in set(drawing_data.loc[(drawing_data[args.compare[0]] == args.compare[1]), "Stage"])) and (x in set(drawing_data.loc[(drawing_data[args.compare[0]] == args.compare[2]), "Stage"])), stage_list))
        drawing_palette = list(map(lambda x: step00.stage_color_code[x], drawing_stage_list))

        try:
            stat, p = scipy.stats.kruskal(*[drawing_data.loc[(drawing_data["Stage"] == stage) & (drawing_data[args.compare[0]] == clinical), watching] for stage, clinical in itertools.product(drawing_stage_list, args.compare[1:])])
        except ValueError:
            stat, p = 0.0, 1.0

        seaborn.violinplot(data=drawing_data, x="Stage", y=watching, order=drawing_stage_list, hue=args.compare[0], hue_order=args.compare[1:], inner="box", palette=drawing_palette, cut=1, ax=ax)
        statannotations.Annotator.Annotator(ax, [((stage, args.compare[1]), (stage, args.compare[2])) for stage in drawing_stage_list], data=drawing_data, x="Stage", y=watching, order=drawing_stage_list, hue=args.compare[0], hue_order=args.compare[1:]).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

        matplotlib.pyplot.title(f"{gene}: Kruskal-Wallis p={p:.3f}")
        matplotlib.pyplot.xlabel("")
        matplotlib.pyplot.ylabel("Ratio")

        matplotlib.pyplot.tight_layout()
        figures.append("{0}.pdf".format(gene))
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
