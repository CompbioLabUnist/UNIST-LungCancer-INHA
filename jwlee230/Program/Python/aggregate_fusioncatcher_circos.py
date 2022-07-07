"""
aggregate_fusioncatcher_circos.py: Aggregate FusionCatcher results as circos plot
"""
import argparse
import multiprocessing
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import tqdm
import step00

size_data = pandas.DataFrame()


def read_maf(filename: str) -> pandas.DataFrame:
    data = pandas.read_csv(filename, sep="\t", low_memory=False)
    data["sample"] = step00.get_id(filename)
    data["stage"] = step00.get_long_sample_type(filename)
    return data


def location_to_angle(chromosome: str, location: int) -> float:
    full_length = sum(size_data["length"])

    current_length = sum(size_data.iloc[:step00.chromosome_full_list.index(chromosome)]["length"])
    current_length += location

    return current_length / full_length * (2 * numpy.pi)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="FusionCatcher output .tsv files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("size", help="Chromosome size", type=str)
    parser.add_argument("output", help="Output file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    size_data = pandas.read_csv(args.size, sep="\t", header=None, names=["chromosome", "length"]).set_index(keys="chromosome", verify_integrity=True)
    size_data = size_data.loc[step00.chromosome_full_list, :]
    print(size_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(sorted(patients))

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))
    args.input.sort(key=step00.sorting_by_type)

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
    input_data["gene_1_chrom"] = list(map(lambda x: "chr" + x.split(":")[0], input_data["Fusion_point_for_gene_1(5end_fusion_partner)"]))
    input_data["gene_1_loc"] = list(map(lambda x: int(x.split(":")[1]), input_data["Fusion_point_for_gene_1(5end_fusion_partner)"]))
    input_data["gene_2_chrom"] = list(map(lambda x: "chr" + x.split(":")[0], input_data["Fusion_point_for_gene_2(3end_fusion_partner)"]))
    input_data["gene_2_loc"] = list(map(lambda x: int(x.split(":")[1]), input_data["Fusion_point_for_gene_2(3end_fusion_partner)"]))
    with multiprocessing.Pool(args.cpus) as pool:
        input_data["gene_1_angle"] = pool.starmap(location_to_angle, input_data[["gene_1_chrom", "gene_1_loc"]].itertuples(index=False, name=None))
        input_data["gene_2_angle"] = pool.starmap(location_to_angle, input_data[["gene_2_chrom", "gene_2_loc"]].itertuples(index=False, name=None))
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24), subplot_kw={"projection": "polar"})

    step = 0.01
    legend_set = set()

    for index, row in tqdm.tqdm(input_data.iterrows()):
        if row["gene_1_angle"] < row["gene_2_angle"]:
            x = numpy.linspace(row["gene_1_angle"], row["gene_2_angle"], step00.small)
        else:
            x = numpy.linspace(row["gene_2_angle"], row["gene_1_angle"], step00.small)

        threshold = abs(row["gene_1_angle"] - row["gene_2_angle"])
        if threshold > numpy.pi / 2:
            threshold -= numpy.pi / 2

        y = 1 - numpy.sin(numpy.linspace(0, numpy.pi, step00.small)) * 0.9 * (threshold * 2 / numpy.pi)

        if row["stage"] not in legend_set:
            matplotlib.pyplot.plot(x, y, color=step00.stage_color_code[row["stage"]], linestyle=step00.stage_linestyle[row["stage"]], linewidth=3, alpha=0.6, label=row["stage"])
            legend_set.add(row["stage"])
        else:
            matplotlib.pyplot.plot(x, y, color=step00.stage_color_code[row["stage"]], linestyle=step00.stage_linestyle[row["stage"]], linewidth=3, alpha=0.6)

    for chromosome, color in tqdm.tqdm(list(zip(step00.chromosome_full_list, matplotlib.colors.XKCD_COLORS))):
        x = numpy.arange(location_to_angle(chromosome, 0), location_to_angle(chromosome, size_data.loc[chromosome, "length"]), step)
        matplotlib.pyplot.plot(x, [1] * len(x), color=color, linewidth=10)
        matplotlib.pyplot.text(location_to_angle(chromosome, size_data.loc[chromosome, "length"] / 2), 1, chromosome, color="k", horizontalalignment="center", verticalalignment="center", fontsize="x-small", bbox={"color": "white", "alpha": 0.5})

    ax.set_theta_direction(-1)
    ax.set_theta_zero_location("N")
    ax.set_rmax(1.1)

    matplotlib.pyplot.grid(False)
    matplotlib.pyplot.xticks([])
    matplotlib.pyplot.yticks([])
    matplotlib.pyplot.legend(loc="upper left")
    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
