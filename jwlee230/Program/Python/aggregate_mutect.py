"""
aggregate_mutect.py: aggregate mutect MAF files
"""
import argparse
import collections
import multiprocessing
from comut import comut
import matplotlib
import pandas
import step00


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("output", help="Output file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--gene", help="Gene number to draw", type=int, default=30)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif args.gene < 1:
        raise ValueError("GENE must be positive!!")

    args.input.sort()
    matplotlib.rcParams.update({"font.family": "serif"})

    my_comut = comut.CoMut()
    my_comut.samples = sorted(list(map(lambda x: x.split("/")[-1].split(".")[0], args.input)), key=step00.sorting)

    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)

    mutect_data = mutect_data.loc[(mutect_data["Variant_Classification"].isin(step00.mutations_list))]
    mutect_data["Tumor_Sample_Barcode"] = list(map(lambda x: x.split(".")[0], mutect_data["Tumor_Sample_Barcode"]))
    mutect_data["Variant_Classification"] = list(map(lambda x: step00.mutations_dict[x], mutect_data["Variant_Classification"]))
    print(mutect_data)

    counter: collections.Counter = collections.Counter(mutect_data["Hugo_Symbol"])
    mutation_data = pandas.DataFrame(counter.most_common(args.gene), columns=["Gene", "Count"])
    print(mutation_data)

    patient_data = pandas.DataFrame()
    patient_data["Tumor_Sample_Barcode"] = my_comut.samples
    patient_data["Patient"] = list(map(lambda x: hash(step00.get_patient(x)), my_comut.samples))
    patient_data["Collection_Type_category"] = "Collection type"
    patient_data["Collection_Type_value"] = list(map(step00.get_long_sample_type, my_comut.samples))
    patient_data["Mutation_Count"] = list(map(lambda x: mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == x)].shape[0], my_comut.samples))
    print(patient_data)

    my_comut.add_sample_indicators(patient_data[["Tumor_Sample_Barcode", "Patient"]].set_axis(labels=step00.sample_columns, axis="columns"), name="Same patient")
    my_comut.add_categorical_data(patient_data[["Tumor_Sample_Barcode", "Collection_Type_category", "Collection_Type_value"]].set_axis(labels=step00.categorical_columns, axis="columns"), name="Collection type", value_order=step00.long_sample_type_list, mapping=step00.collection_mapping)
    my_comut.add_categorical_data(mutect_data[["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification"]].set_axis(labels=step00.categorical_columns, axis="columns"), name="Mutation type", category_order=mutation_data["Gene"], mapping=step00.mutation_mapping, value_order=["Missense"], priority=["Missense"])
    my_comut.add_bar_data(patient_data[["Tumor_Sample_Barcode", "Mutation_Count"]].set_axis(labels=step00.sample_columns, axis="columns"), name="Mutation count", ylabel="Counts", mapping={"group": "purple"})
    my_comut.add_side_bar_data(mutation_data.set_axis(labels=step00.bar_columns, axis="columns"), name="Mutation count", xlabel="Counts", paired_name="Mutation type")

    my_comut.plot_comut(x_padding=0.04, y_padding=0.04, tri_padding=0.03, figsize=(32, 18))
    my_comut.add_unified_legend()
    my_comut.figure.savefig(args.output)
