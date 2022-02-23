"""
draw_venn_plot_gistic.py: draw venn diagram plot with gistic peak
"""
import argparse
import multiprocessing
import re
import typing
import matplotlib
import matplotlib.pyplot
import pandas
import venn
import tqdm
import step00

band_data = pandas.DataFrame()
cgc_data = pandas.DataFrame()


def generate_logics(n_sets):
    for i in range(1, 2 ** n_sets):
        yield bin(i)[2:].zfill(n_sets)


def sorting_cytoband(s):
    s_split = list(map(int, re.split(r"p|q|\.", s)))

    if "p" in s:
        s_split.insert(1, "p")
    elif "q" in s:
        s_split.insert(1, "q")
    else:
        raise Exception("Something went wrong!!")

    return s_split


def get_chromosome(location: str) -> str:
    return "chr" + location.split(":")[0]


def get_start(location: str) -> int:
    return int(location.replace("-", ":").split(":")[1])


def get_end(location: str) -> int:
    return int(location.replace("-", ":").split(":")[2])


def query_band(chromosome: str, start: int, end: int) -> typing.List[str]:
    answer = list()

    answer += list(band_data.loc[(band_data["chrom"] == chromosome) & (band_data["chrom_start"] <= start) & (start <= band_data["chrom_end"]), "chrom-arm"])
    answer += list(band_data.loc[(band_data["chrom"] == chromosome) & (start <= band_data["chrom_start"]) & (band_data["chrom_end"] <= end), "chrom-arm"])
    answer += list(band_data.loc[(band_data["chrom"] == chromosome) & (band_data["chrom_start"] <= end) & (end <= band_data["chrom_end"]), "chrom-arm"])

    return sorted(set(answer))


def query_gene(cytoband: str) -> str:
    d = list(cgc_data.loc[(cgc_data["Chromosome-Arm"] == cytoband), "Gene Symbol"])
    if d:
        if len(d) < 10:
            return ",".join(d)
        else:
            return "{0},...({1})".format(d[0], len(d))
    else:
        return ""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input Gistic all_lesions.conf_99.txt file(s)", type=str, nargs="+")
    parser.add_argument("cgc", help="CGC CSV files", type=str)
    parser.add_argument("band", help="Chromosome band txt file", type=str)
    parser.add_argument("output", help="Output file basename", type=str)
    parser.add_argument("--annotation", help="Annotation for venn diagram", type=str, nargs="+", required=True)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.05)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--amplification", help="Draw Amplification pathway", action="store_true", default=False)
    group.add_argument("--deletion", help="Draw Deletion pathway", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith("/all_lesions.conf_99.txt"), args.input)):
        raise ValueError("Input is not valid!!")
    elif len(args.input) != len(args.annotation):
        raise ValueError("Annotation must be one-to-one upon DEG!!")
    elif not args.cgc.endswith(".csv"):
        raise ValueError("One must end with .CSV!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.p < 1):
        raise ValueError("P-value must be between 0 and 1!!")

    band_data = step00.get_band_data(args.band)
    band_data["chrom-arm"] = list(map(lambda x: "{0}{1}".format(x[0][3:], x[1]), band_data[["chrom", "name"]].itertuples(index=False, name=None)))
    chromosome_list = sorted(set(band_data["chrom-arm"]))
    print(band_data)

    cgc_data = pandas.read_csv(args.cgc)
    cgc_data = cgc_data.loc[~(cgc_data["Genome Location"].str.contains(":-"))]
    with multiprocessing.Pool(args.cpus) as pool:
        cgc_data["Chromosome"] = pool.map(get_chromosome, cgc_data["Genome Location"])
        cgc_data["Start"] = pool.map(get_start, cgc_data["Genome Location"])
        cgc_data["End"] = pool.map(get_end, cgc_data["Genome Location"])
        cgc_data["Chromosome-Arm"] = pool.starmap(query_band, cgc_data[["Chromosome", "Start", "End"]].itertuples(index=False, name=None))
    cgc_data = cgc_data.explode("Chromosome-Arm", ignore_index=True)
    print(cgc_data)

    input_data = dict()
    for annotation, input_file in tqdm.tqdm(list(zip(args.annotation, args.input))):
        data = pandas.read_csv(input_file, sep="\t")
        data = data.loc[(data["q values"] < args.p) & (data["Residual q values after removing segments shared with higher peaks"] < args.p), :]
        data["Descriptor"] = list(map(lambda x: x.strip(), data["Descriptor"]))
        if args.amplification:
            input_data[annotation] = set(data.loc[(data["Unique Name"].str.contains("Amplification Peak")) & ~(data["Unique Name"].str.contains("CN")), "Descriptor"])
        elif args.deletion:
            input_data[annotation] = set(data.loc[(data["Unique Name"].str.contains("Deletion Peak")) & ~(data["Unique Name"].str.contains("CN")), "Descriptor"])
        else:
            raise Exception("Something went wrong!!")
        input_data[annotation] = set(filter(lambda x: (not x.startswith("X")) and (not x.startswith("Y")), input_data[annotation]))
    print(input_data)

    every_cytoband = sorted(set.union(*list(input_data.values())), key=sorting_cytoband)
    print(every_cytoband)

    if every_cytoband:
        output_data = pandas.DataFrame(data=[["" for x in args.annotation] for y in every_cytoband], index=every_cytoband, columns=args.annotation, dtype=str)
        for annotation in tqdm.tqdm(args.annotation):
            output_data.loc[input_data[annotation], annotation] = "*"
    else:
        output_data = pandas.DataFrame(data=[["" for x in args.annotation]], index=[""], columns=args.annotation, dtype=str)

    with multiprocessing.Pool(args.cpus) as pool:
        output_data["Gene"] = pool.map(query_gene, list(output_data.index))

    print(output_data)
    output_data.to_latex(args.output + ".tex", column_format="l" + "c" * len(args.annotation) + "r")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    try:
        venn.venn(input_data, ax=ax, fmt=step00.venn_format, fontsize=step00.matplotlib_parameters["legend.fontsize"], legend_loc="upper left")
    except ZeroDivisionError:
        matplotlib.pyplot.text(0.5, 0.5, "Nothing to show...", fontsize=step00.matplotlib_parameters["axes.titlesize"], color="k", horizontalalignment="center", verticalalignment="center")
        matplotlib.pyplot.xticks([])
        matplotlib.pyplot.yticks([])
    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output + ".pdf")
    matplotlib.pyplot.close(fig)

    datasets = list(input_data.values())
    dataset_union = set(every_cytoband)
    for logic in generate_logics(len(args.annotation)):
        included_sets = [datasets[i] for i in range(len(args.annotation)) if logic[i] == "1"]
        excluded_sets = [datasets[i] for i in range(len(args.annotation)) if logic[i] == "0"]

        print("+".join(list(map(lambda y: y[0], list(filter(lambda x: x[1] == "1", zip(args.annotation, list(logic))))))))
        print(sorted((dataset_union & set.intersection(*included_sets)) - set.union(set(), *excluded_sets), key=sorting_cytoband))
        print()
