"""
aggregate_mutect.py: aggregate mutect MAF files
"""
import argparse
import collections
import multiprocessing
from comut import comut
from comut import fileparsers
import matplotlib
import matplotlib.pyplot
import pandas
import step00


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("output", help="Output file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    args.input.sort()

    mut_mapping = {"Missense": "green", "Nonsense": "deeppink", "In frame indel": {"facecolor": "blue", "hatch": "xxx"}, "Frameshift indel": "#FFD700", "Splice site": "darkviolet", "LOH": {"facecolor": "none", "edgecolor": "black", "linewidth": 3}, "Absent": {"facecolor": "white", "alpha": 0.2}}

    example_comut = comut.CoMut()
    example_comut.samples = list(map(lambda x: x.split("/")[-1].split(".")[0], args.input))

    with multiprocessing.Pool(args.cpus) as pool:
        mutation_df = fileparsers.parse_maf(pandas.concat(pool.map(read_maf, args.input[:10]), ignore_index=True, copy=False))
    mutation_df["sample"] = list(map(lambda x: x.split(".")[0], mutation_df["sample"]))
    print(mutation_df)

    category_counter: collections.Counter = collections.Counter(mutation_df["category"])
    important_mutations = list(map(lambda x: x[0], category_counter.most_common(10)))

    example_comut.add_categorical_data(mutation_df, name="Mutation type", category_order=important_mutations, mapping=mut_mapping)
    example_comut.plot_comut(x_padding=0.04, y_padding=0.04, tri_padding=0.03, figsize=(32, 18))
    example_comut.figure.savefig(args.output)
