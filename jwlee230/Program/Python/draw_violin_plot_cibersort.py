"""
draw_violin_cibersort.py: draw violin plot from CibersortX result
"""
import argparse
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import statannot
import step00

cibersort_data = pandas.DataFrame()


def run(cell: str) -> str:
    print("Running:", cell)
    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
    tmp = set(cibersort_data["Stage"])

    order = list(filter(lambda x: x in tmp, step00.long_sample_type_list))

    seaborn.violinplot(data=cibersort_data, x="Stage", y=cell, order=order)
    statannot.add_stat_annotation(ax, data=cibersort_data, x="Stage", y=cell, order=order, test="Mann-Whitney", box_pairs=itertools.combinations(order, 2), text_format="star", loc="outside", verbose=0)

    matplotlib.pyplot.title(cell)
    matplotlib.pyplot.ylabel("Proportion")
    matplotlib.pyplot.tight_layout()

    fig_name = cell.replace(" ", "").replace("(", "").replace(")", "").replace("/", "_") + ".pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("cibersort", help="CIBERSORT result TSV file (not necessarily TSV)", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if not args.output.endswith(".tar"):
        raise ValueError("Output must end with .tar!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    cibersort_data = pandas.read_csv(args.cibersort, sep="\t", index_col="Mixture")
    cibersort_data["Stage"] = list(map(step00.get_long_sample_type, list(cibersort_data.index)))
    print(cibersort_data)

    with multiprocessing.Pool(args.cpus) as pool:
        tar_files = sorted(pool.map(run, list(cibersort_data.columns)[:-1]))

    with tarfile.open(args.output, "w") as tar:
        for f in tar_files:
            tar.add(f, arcname=f)
