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


def run(cell: str, ADC: bool = False, SQC: bool = False) -> str:
    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    if ADC:
        order = step00.simple_stage_list
    elif SQC:
        order = step00.simple_stage_list
    else:
        raise Exception("Something went wrong!!")

    seaborn.violinplot(data=cibersort_data, x="Stage", y=cell, order=order)
    statannot.add_stat_annotation(ax, data=cibersort_data, x="Stage", y=cell, order=order, test="t-test_ind", box_pairs=itertools.combinations(order, 2), text_format="star", loc="inside", verbose=0)

    matplotlib.pyplot.ylabel("Proportion")
    if ADC:
        matplotlib.pyplot.title("{0} in ADC".format(cell.title()))
    elif SQC:
        matplotlib.pyplot.title("{0} in SQC".format(cell.title()))
    else:
        raise Exception("Something went wrong!!")
    ax.figure.tight_layout()

    fig_name = cell.title().replace(" ", "").replace("(", "").replace(")", "") + ".pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("count", help="Count TSV file", type=str)
    parser.add_argument("cibersort", help="CIBERSORT result TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ADC", help="Draw ADC pathway", action="store_true", default=False)
    group.add_argument("--SQC", help="Draw SQC pathway", action="store_true", default=False)

    args = parser.parse_args()

    if not args.count.endswith(".tsv"):
        raise ValueError("Count must end with .TSV!!")
    elif not args.cibersort.endswith(".tsv"):
        raise ValueError("CIBERSORT must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .tar!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    count_data = pandas.read_csv(args.count, sep="\t", index_col="gene_name")
    print(count_data)

    cibersort_data = pandas.read_csv(args.cibersort, sep="\t", index_col="Column")
    cibersort_data["Patient_ID"] = list(count_data.columns)
    cibersort_data.set_index("Patient_ID", inplace=True)
    cibersort_data.drop(columns=list(cibersort_data.columns)[-4:], inplace=True)
    cell_list = list(cibersort_data.columns)
    cibersort_data["Stage"] = list(map(step00.get_simple_sample_type, list(cibersort_data.index)))
    print(cibersort_data)

    with multiprocessing.Pool(args.cpus) as pool:
        tar_files = pool.starmap(run, [(cell, args.ADC, args.SQC) for cell in cell_list])

    with tarfile.open(args.output, "w") as tar:
        for f in tar_files:
            tar.add(f, arcname=f)
