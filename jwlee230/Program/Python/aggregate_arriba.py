"""
aggregate_arriba.py: Aggregate the Arriba resuluts as CoMut plot
"""
import argparse
import collections
import multiprocessing
from comut import comut
import matplotlib
import pandas
import step00


def read_maf(filename: str) -> pandas.DataFrame:
    data = pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)
    data = data.loc[(data["confidence"] == "high")]
    data["sample"] = filename.split("/")[-1].split(".")[0]
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Arriba output .tsv files", type=str, nargs="+")
    parser.add_argument("output", help="Output file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--gene", help="Gene number to draw", type=int, default=50)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif args.gene < 1:
        raise ValueError("GENE must be positive!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    my_comut = comut.CoMut()
    my_comut.samples = sorted(list(map(lambda x: x.split("/")[-1].split(".")[0], args.input)), key=step00.sorting)

    with multiprocessing.Pool(args.cpus) as pool:
        arriba_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
    print(arriba_data)

    gene_counter = collections.Counter(list(arriba_data["gene1"]) + list(arriba_data["gene2"]))
    common_genes = set(map(lambda x: x[0], gene_counter.most_common(args.gene)))
    print(gene_counter.most_common(args.gene))
    print(common_genes)

    indicator_data = pandas.DataFrame()
    indicator_data["sample"] = my_comut.samples
    indicator_data["group"] = list(map(lambda x: hash(step00.get_patient(x)), indicator_data["sample"]))
    print(indicator_data)
    my_comut.add_sample_indicators(indicator_data, name="Same patient")

    categorical_data = pandas.DataFrame()
    categorical_data["sample"] = arriba_data["sample"]
    categorical_data["value"] = list(map(lambda x: x.split("/")[0].title(), arriba_data["type"]))
    categorical_data["category"] = list(map(lambda x: "{0}({2}) + {1}({3})".format(x[0], x[1], x[2], x[3]), arriba_data[["gene1", "gene2", "breakpoint1", "breakpoint2"]].itertuples(index=False, name=None)))
    # categorical_data = categorical_data.loc[(arriba_data["gene1"].isin(common_genes)) & (arriba_data["gene2"].isin(common_genes))]
    print(categorical_data)
    my_comut.add_categorical_data(categorical_data, name="Gene fusion type")

    bar_data = pandas.DataFrame()
    bar_data["sample"] = my_comut.samples
    bar_data["Number of fusion"] = list(map(lambda x: len(arriba_data[(arriba_data["sample"] == x)]), bar_data["sample"]))
    print(bar_data)
    my_comut.add_bar_data(bar_data, name="Number of fusion", ylabel="Counts")

    my_comut.plot_comut(x_padding=0.04, y_padding=0.04, tri_padding=0.03, figsize=(len(my_comut.samples), len(categorical_data)))
    my_comut.add_unified_legend()
    my_comut.figure.savefig(args.output, bbox_inches="tight")
