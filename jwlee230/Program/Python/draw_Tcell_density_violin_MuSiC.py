"""
draw_Tcell_density_violin_MuSiC.py: draw T-cell density violin plot for MuSiC
"""
import argparse
import itertools
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import statannotations.Annotator
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="BisqueRNA result TSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .tsv!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .pdf!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)
    input_data = input_data.reindex(index=sorted(list(input_data.index), key=step00.sorting_by_type)).T
    input_data = input_data.reindex(index=sorted(list(input_data.index), key=lambda x: sum(input_data.loc[x, :]), reverse=True))
    print(input_data)
