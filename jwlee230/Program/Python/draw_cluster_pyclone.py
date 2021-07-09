"""
draw_cluster_pyclone.py: draw lineplot upon cluster results by PyClone
"""
import argparse
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("cluster", help="Cluster TSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    args = parser.parse_args()

    if not args.cluster.endswith(".tsv"):
        raise ValueError("Cluster must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    cluster_data = pandas.read_csv(args.cluster, sep="\t")
    print(cluster_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))
    seaborn.lineplot(data=cluster_data, x="sample_id", y="mean", hue="cluster_id", markers=True)

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
