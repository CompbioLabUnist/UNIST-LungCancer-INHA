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
import tqdm
import step00

sample_dict: typing.Dict[str, typing.Dict[str, typing.List[str]]] = dict()
gene_set = set()


def draw_plot(first_sample: str, second_sample: str) -> typing.List[str]:
    first_name = step00.get_id(first_sample)
    second_name = step00.get_id(second_sample)

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
    merged_data["mutation"] = list(map(lambda x: x[7] in step00.nonsynonymous_mutations, list(merged_data.index)))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    texts = list()

    data = merged_data.loc[~(merged_data["gene_census"]) & ~(merged_data["mutation"])]
    matplotlib.pyplot.scatter(data["first_VAF"], data["second_VAF"], c="tab:gray", marker="o", alpha=0.3, s=12 ** 2, edgecolor="none", label="Synonymous mutations")

    data = merged_data.loc[~(merged_data["gene_census"]) & (merged_data["mutation"])]
    matplotlib.pyplot.scatter(data["first_VAF"], data["second_VAF"], c="black", marker="*", alpha=0.3, s=15 ** 2, edgecolor="none", label="Functional mutations")

    data = merged_data.loc[(merged_data["gene_census"]) & (merged_data["mutation"])]
    matplotlib.pyplot.scatter(data["first_VAF"], data["second_VAF"], c="tab:red", marker="*", alpha=1.0, s=20 ** 2, edgecolor="none", label="Cancer genes")

    for index, d in data.iterrows():
        if (d["first_VAF"] > 0.6) and (d["second_VAF"] == 0.0):
            texts.append(matplotlib.pyplot.text(d["first_VAF"], d["second_VAF"], "{0}: {1}".format(index[6], index[8]), fontsize="xx-small", bbox={"color": "white", "alpha": 0.5}))
        elif (d["first_VAF"] == 0.0) and (d["second_VAF"] > 0.6):
            texts.append(matplotlib.pyplot.text(d["first_VAF"], d["second_VAF"], "{0}: {1}".format(index[6], index[8]), fontsize="xx-small", bbox={"color": "white", "alpha": 0.5}))
        elif (d["first_VAF"] > 0.0) and (d["second_VAF"] > 0.0):
            texts.append(matplotlib.pyplot.text(d["first_VAF"], d["second_VAF"], "{0}: {1}".format(index[6], index[8]), fontsize="small", bbox={"color": "white", "alpha": 0.5}))

    matplotlib.pyplot.axline((0, 0), (1, 1), linestyle="--", color="black", alpha=0.3)
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.xlim(-0.1, 1.1)
    matplotlib.pyplot.ylim(-0.1, 1.1)
    matplotlib.pyplot.xlabel(f"VAF of {first_name} ({step00.get_long_sample_type(first_name)})")
    matplotlib.pyplot.ylabel(f"VAF of {second_name} ({step00.get_long_sample_type(second_name)})")
    matplotlib.pyplot.title(f"{first_name} vs. {second_name}")
    matplotlib.pyplot.legend(loc="upper right")
    matplotlib.pyplot.tight_layout()
    adjust_text(texts, arrowprops={"arrowstyle": "-", "color": "k", "linewidth": 0.5}, ax=ax, lim=step00.big)

    figure_name = f"{first_name}+{second_name}.pdf"
    fig.savefig(figure_name)
    matplotlib.pyplot.close(fig)

    tsv_name = f"{first_name}+{second_name}.tsv"
    merged_data.to_csv(tsv_name, sep="\t", header=True, index=True)

    return [figure_name, tsv_name]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data data CSV file", type=str)
    parser.add_argument("gene", help="Cancer gene census CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_histology = parser.add_mutually_exclusive_group(required=True)
    group_histology.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_histology.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

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

    patients = list(map(lambda x: step00.get_patient(x.split("/")[-1].split(".")[0]), args.input))

    if args.SQC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
        patients = list(filter(lambda x: step00.get_patient(x) in histology, patients))
    elif args.ADC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
        patients = list(filter(lambda x: step00.get_patient(x) in histology, patients))
    else:
        raise Exception("Something went wrong!!")
    print(patients)
    args.input = list(filter(lambda x: step00.get_patient(step00.get_id(x)), args.input))

    gene_data = pandas.read_csv(args.gene)
    gene_set = set(gene_data["Gene Symbol"])

    for input_file in tqdm.tqdm(args.input):
        cleared_input_file = input_file.split("/")[-1].split(".")[0]
        patient = step00.get_patient(cleared_input_file)
        stage = step00.get_long_sample_type(cleared_input_file)

        if patient not in sample_dict:
            sample_dict[patient] = dict()

        if stage not in sample_dict[patient]:
            sample_dict[patient][stage] = list()

        sample_dict[patient][stage].append(input_file)
    pprint.pprint(sample_dict)

    compare_list = list()
    for patient in tqdm.tqdm(sample_dict):
        for first_type, second_type in itertools.combinations(step00.long_sample_type_list, 2):
            if (first_type not in sample_dict[patient]) or (second_type not in sample_dict[patient]):
                continue
            compare_list += list(itertools.product(sample_dict[patient][first_type], sample_dict[patient][second_type]))
    print(len(compare_list))
    pprint.pprint(compare_list)

    with multiprocessing.Pool(args.cpus) as pool:
        files = sorted(itertools.chain.from_iterable(pool.starmap(draw_plot, compare_list)))

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(files):
            tar.add(f, arcname=f)
