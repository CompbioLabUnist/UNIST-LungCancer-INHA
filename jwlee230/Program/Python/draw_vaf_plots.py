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


def draw_plot(first_sample: str, second_sample: str) -> typing.List[str]:
    first_name = first_sample.split("/")[-1].split(".")[0]
    second_name = second_sample.split("/")[-1].split(".")[0]

    first_data = pandas.read_csv(first_sample, sep="\t", comment="#", low_memory=False)
    second_data = pandas.read_csv(second_sample, sep="\t", comment="#", low_memory=False)

    first_data["first_VAF"] = first_data["t_alt_count"] / first_data["t_depth"]
    second_data["second_VAF"] = second_data["t_alt_count"] / second_data["t_depth"]

    first_data = first_data.loc[(first_data["Chromosome"].isin(step00.chromosome_list)), :]
    second_data = second_data.loc[(second_data["Chromosome"].isin(step00.chromosome_list)), :]

    first_data.set_index(keys=["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Hugo_Symbol", "Variant_Classification", "HGVSp_Short"], inplace=True, verify_integrity=True)
    second_data.set_index(keys=["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Hugo_Symbol", "Variant_Classification", "HGVSp_Short"], inplace=True, verify_integrity=True)

    first_data = first_data[["first_VAF"]]
    second_data = second_data[["second_VAF"]]

    merged_data = pandas.concat(objs=[first_data, second_data], axis="columns", join="outer", sort=True, copy=False).fillna(value=0.0)
    merged_data["gene_census"] = list(map(lambda x: x[6] in gene_set, list(merged_data.index)))
    merged_data["mutation"] = list(map(lambda x: x[7] in step00.mutations_list, list(merged_data.index)))

    print("{0} vs {1}: {2}".format(first_name, second_name, merged_data.shape))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    texts = []

    data = merged_data.loc[~(merged_data["gene_census"]) & ~(merged_data["mutation"])]
    matplotlib.pyplot.scatter(data["first_VAF"], data["second_VAF"], c="tab:gray", marker="o", alpha=0.3, s=12 ** 2, edgecolor="none", label="Synonymous mutations")

    data = merged_data.loc[~(merged_data["gene_census"]) & (merged_data["mutation"])]
    matplotlib.pyplot.scatter(data["first_VAF"], data["second_VAF"], c="black", marker="*", alpha=0.3, s=15 ** 2, edgecolor="none", label="Functional mutations")

    data = merged_data.loc[(merged_data["gene_census"]) & (merged_data["mutation"])]
    matplotlib.pyplot.scatter(data["first_VAF"], data["second_VAF"], c="tab:red", marker="*", alpha=1.0, s=20 ** 2, edgecolor="none", label="Cancer genes")
    for index, d in data.iterrows():
        if (d["first_VAF"] > 0.6) and (d["second_VAF"] == 0.0):
            texts.append(matplotlib.pyplot.text(d["first_VAF"], d["second_VAF"], "{0}: {1}".format(index[6], index[8]), fontsize="xx-small"))
        elif (d["first_VAF"] == 0.0) and (d["second_VAF"] > 0.6):
            texts.append(matplotlib.pyplot.text(d["first_VAF"], d["second_VAF"], "{0}: {1}".format(index[6], index[8]), fontsize="xx-small"))
        elif (d["first_VAF"] > 0.0) and (d["second_VAF"] > 0.0):
            texts.append(matplotlib.pyplot.text(d["first_VAF"], d["second_VAF"], "{0}: {1}".format(index[6], index[8]), fontsize="small"))

    matplotlib.pyplot.axline((0, 0), (1, 1), linestyle="--", color="black", alpha=0.3)
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.xlim(-0.1, 1.1)
    matplotlib.pyplot.ylim(-0.1, 1.1)
    matplotlib.pyplot.xlabel("VAF of {0} ({1})".format(first_name, step00.get_long_sample_type(first_name)))
    matplotlib.pyplot.ylabel("VAF of {0} ({1})".format(second_name, step00.get_long_sample_type(second_name)))
    matplotlib.pyplot.title("{0} vs. {1}".format(first_name, second_name))
    matplotlib.pyplot.legend(loc="upper right")
    adjust_text(texts, arrowprops={"arrowstyle": "-", "color": "k", "linewidth": 0.5}, ax=ax, lim=10 ** 6)

    figure_name = "{0}+{1}.pdf".format(first_name, second_name)
    fig.savefig(figure_name)
    matplotlib.pyplot.close(fig)

    tsv_name = "{0}+{1}.tsv".format(first_name, second_name)
    merged_data.to_csv(tsv_name, sep="\t", header=True, index=True)

    return [figure_name, tsv_name]


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
        files = itertools.chain.from_iterable(pool.starmap(draw_plot, compare_list))

    with tarfile.open(args.output, "w") as tar:
        for f in sorted(files):
            print("Compressing:", f)
            tar.add(f, arcname=f)
