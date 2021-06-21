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

    first_data = pandas.read_csv(first_sample, sep="\t", comment="#", low_memory=False, index_col=["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Hugo_Symbol", "Variant_Classification", "HGVSp_Short"])
    second_data = pandas.read_csv(second_sample, sep="\t", comment="#", low_memory=False, index_col=["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Hugo_Symbol", "Variant_Classification", "HGVSp_Short"])

    first_data["first_VAF"] = first_data["t_alt_count"] / first_data["t_depth"]
    second_data["second_VAF"] = second_data["t_alt_count"] / second_data["t_depth"]

    first_data = first_data[["first_VAF"]]
    second_data = second_data[["second_VAF"]]

    merged_data = pandas.concat(objs=[first_data, second_data], axis="columns", join="outer", sort=True, copy=False).fillna(value=0.0)
    merged_data["gene_census"] = list(map(lambda x: x[5] in gene_set, list(merged_data.index)))
    merged_data["mutation"] = list(map(lambda x: x[6] in step00.mutations_list, list(merged_data.index)))

    print("{0} vs {1}: {2}".format(first_name, second_name, merged_data.shape))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    texts = []

    data = merged_data.loc[~(merged_data["gene_census"]) & ~(merged_data["mutation"])]
    matplotlib.pyplot.scatter(data["first_VAF"], data["second_VAF"], c="tab:gray", marker="o", alpha=0.3, s=12 ** 2, edgecolor="none", label="Synonymous Mutations")

    data = merged_data.loc[~(merged_data["gene_census"]) & (merged_data["mutation"])]
    matplotlib.pyplot.scatter(data["first_VAF"], data["second_VAF"], c="black", marker="*", alpha=0.3, s=12 ** 2, edgecolor="none", label="SNPs")

    data = merged_data.loc[(merged_data["gene_census"]) & (merged_data["mutation"])]
    matplotlib.pyplot.scatter(data["first_VAF"], data["second_VAF"], c="tab:red", marker="*", alpha=1.0, s=20 ** 2, edgecolor="none", label="Meaningful genes")
    for index, d in data.iterrows():
        if (d["first_VAF"] == 0) or (d["second_VAF"] == 0):
            continue
        texts.append(matplotlib.pyplot.text(d["first_VAF"], d["second_VAF"], "{0}: {1}".format(index[5], index[7]), fontsize="small"))

    matplotlib.pyplot.axline((0, 0), (1, 1), linestyle="--", color="black", alpha=0.3)
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.xlim(-0.1, 1.1)
    matplotlib.pyplot.ylim(-0.1, 1.1)
    matplotlib.pyplot.xlabel("VAF of {0} ({1})".format(first_name, step00.get_long_sample_type(first_name)))
    matplotlib.pyplot.ylabel("VAF of {0} ({1})".format(second_name, step00.get_long_sample_type(second_name)))
    matplotlib.pyplot.title("{0} vs. {1}".format(first_name, second_name))
    matplotlib.pyplot.legend(loc="upper right")
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
