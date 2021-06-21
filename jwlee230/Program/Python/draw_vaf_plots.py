"""
draw_vaf_plots.py: Draw VAF plots upon Mutect2 MAF results
"""
import argparse
import itertools
import multiprocessing
import tarfile
import typing
import pprint
from adjustText import adjust_text
import matplotlib
import matplotlib.pyplot
import pandas
import step00

sample_dict: typing.Dict[str, typing.Dict[str, typing.List[str]]] = dict()
gene_set = set()


def draw_plot(first_sample: str, second_sample: str) -> str:
    first_name = first_sample.split("/")[-1].split(".")[0]
    second_name = second_sample.split("/")[-1].split(".")[0]

    first_data = pandas.read_csv(first_sample, sep="\t", comment="#", low_memory=False)
    second_data = pandas.read_csv(second_sample, sep="\t", comment="#", low_memory=False)

    first_data["VAF"] = first_data["t_alt_count"] / first_data["t_depth"]
    second_data["VAF"] = second_data["t_alt_count"] / second_data["t_depth"]

    intersected_positions = sorted(set(first_data.loc[:, ["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]].itertuples(index=False, name=None)) | set(second_data.loc[:, ["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]].itertuples(index=False, name=None)))
    print("{0} vs {1}: {2} genes".format(first_name, second_name, len(intersected_positions)))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    texts = []
    for chromosome, start_position, end_position, ref_allele, tumor_seq_allele1, tumor_seq_allele2 in intersected_positions:
        if (d := first_data.loc[(first_data["Chromosome"] == chromosome) & (first_data["Start_Position"] == start_position) & (first_data["End_Position"] == end_position), ["Hugo_Symbol", "Variant_Classification", "VAF", "HGVSp_Short"]]).empty:
            first_row = []
            first_vaf = 0
        else:
            first_row = list(d.to_numpy()[0])
            first_vaf = first_row[2]

        if (d := second_data.loc[(second_data["Chromosome"] == chromosome) & (second_data["Start_Position"] == start_position) & (second_data["End_Position"] == end_position), ["Hugo_Symbol", "Variant_Classification", "VAF", "HGVSp_Short"]]).empty:
            second_row = []
            second_vaf = 0
        else:
            second_row = list(d.to_numpy()[0])
            second_vaf = second_row[2]

        if first_row:
            symbol, variant, mutation = first_row[0], first_row[1], first_row[3]
        else:
            symbol, variant, mutation = second_row[0], second_row[1], second_row[3]

        if (symbol in gene_set) and (variant in step00.mutations_list):
            c = "tab:red"
            marker = "*"
            alpha = 1.0
            s = 20 ** 2
            texts.append(matplotlib.pyplot.text(first_vaf, second_vaf, "{0}: {1}".format(symbol, mutation), fontsize="small"))
        elif (variant in step00.mutations_list):
            c = "tab:gray"
            marker = "*"
            alpha = 0.7
            s = 12 ** 2
        else:
            c = "tab:gray"
            marker = "o"
            alpha = 0.3
            s = 12 ** 2

        matplotlib.pyplot.scatter(first_vaf, second_vaf, c=c, marker=marker, alpha=alpha, s=s, edgecolor="none")

    matplotlib.pyplot.axline((0, 0), (1, 1), linestyle="--", color="black", alpha=0.3)
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.xlim(-0.1, 1.1)
    matplotlib.pyplot.ylim(-0.1, 1.1)
    matplotlib.pyplot.xlabel("VAF of {0} ({1})".format(first_name, step00.get_long_sample_type(first_name)))
    matplotlib.pyplot.ylabel("VAF of {0} ({1})".format(second_name, step00.get_long_sample_type(second_name)))
    matplotlib.pyplot.title("{0} vs. {1}".format(first_name, second_name))
    adjust_text(texts, arrowprops={"arrowstyle": "-", "color": "k", "linewidth": 0.5}, ax=ax, lim=10 ** 5)

    figure_name = "{0}+{1}.pdf".format(first_name, second_name)
    fig.savefig(figure_name)
    matplotlib.pyplot.close(fig)
    return figure_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("gene", help="Cancer gene census CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    gene_data = pandas.read_csv(args.gene)
    gene_set = set(gene_data["Gene Symbol"])

    for input_file in args.input:
        cleared_input_file = input_file.split("/")[-1].split(".")[0]
        patient = step00.get_patient(cleared_input_file)
        stage = step00.get_long_sample_type(cleared_input_file)

        if patient not in sample_dict:
            sample_dict[patient] = dict()

        if stage not in sample_dict[patient]:
            sample_dict[patient][stage] = list()

        sample_dict[patient][stage].append(input_file)
    pprint.pprint(sample_dict)

    compare_list = []
    for patient in sample_dict:
        for first_type, second_type in itertools.combinations(step00.long_sample_type_list, 2):
            if (first_type not in sample_dict[patient]) or (second_type not in sample_dict[patient]):
                continue
            compare_list += list(itertools.product(sample_dict[patient][first_type], sample_dict[patient][second_type]))

    print(len(compare_list))

    with multiprocessing.Pool(args.cpus) as pool:
        figures = pool.starmap(draw_plot, compare_list)

    with tarfile.open(args.output, "w") as tar:
        for f in figures:
            tar.add(f, arcname=f)
