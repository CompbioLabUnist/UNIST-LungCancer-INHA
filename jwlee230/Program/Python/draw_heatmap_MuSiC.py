"""
draw_heatmap_MuSiC.py: draw heatmap plot from MuSiC result
"""
import argparse
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="MuSiC result TSV file ", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .tsv!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    input_data = pandas.read_csv(args.input, sep="\t")
    input_data = input_data.set_index(list(input_data.columns)[0])
    print(input_data)
    input_data = input_data.reindex(index=sorted(list(input_data.index), key=step00.sorting_by_type)).T
    input_data = input_data.reindex(index=sorted(list(input_data.index), key=lambda x: sum(input_data.loc[x, :]), reverse=True))
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(input_data.shape[1], input_data.shape[0]))

    seaborn.heatmap(data=input_data, xticklabels=True, yticklabels=True, cbar=False, square=False, ax=ax, cmap="Reds", vmin=0, vmax=1)

    matplotlib.pyplot.xticks(fontsize="xx-small")
    matplotlib.pyplot.yticks(fontsize="xx-small")
    ax.figure.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
