"""
compare_gistic_peaks.py: compare gistic peaks
"""
import argparse
import typing
import matplotlib
import matplotlib.pyplot
import pandas
import tqdm
import upsetplot
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input Gistic [amp_genes.conf_99.txt|del_genes.conf_99.txt] file(s)", type=str, nargs="+")
    parser.add_argument("cgc", help="CGC CSV files", type=str)
    parser.add_argument("output", help="Output file basename", type=str)
    parser.add_argument("--annotation", help="Annotation for venn diagram", type=str, nargs="+", required=True)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.05)

    args = parser.parse_args()

    if list(filter(lambda x: not (x.endswith("/amp_genes.conf_99.txt") or x.endswith("/del_genes.conf_99.txt")), args.input)):
        raise ValueError("Input is not valid!!")
    elif len(args.input) != len(args.annotation):
        raise ValueError("Annotation must be one-to-one upon DEG!!")
    elif not args.cgc.endswith(".csv"):
        raise ValueError("One must end with .CSV!!")
    elif not (0 < args.p < 1):
        raise ValueError("P-value must be between 0 and 1!!")

    cgc_data = pandas.read_csv(args.cgc)
    cgc_genes = set(cgc_data["Gene Symbol"])
    print(cgc_data)

    input_data: typing.Dict[str, typing.Set[str]] = dict()
    for annotation, input_file in tqdm.tqdm(list(zip(args.annotation, args.input))):
        input_data[annotation] = set()

        data = pandas.read_csv(input_file, sep="\t", keep_default_na=False)
        for column in list(data.columns)[1:-1]:
            if column.startswith("X") or column.startswith("Y"):
                continue
            elif (float(data.loc[0, column]) > args.p) or (float(data.loc[1, column]) > args.p):
                continue

            input_data[annotation] |= set(filter(None, list(map(lambda x: x.strip("[]"), data.loc[3:, column]))))
    print(input_data)

    every_genes = sorted(set.union(*list(input_data.values())))
    index_name = "Genes"
    if set(every_genes) & cgc_genes:
        index_name = "CGC Genes"
        every_genes = sorted(set(every_genes) & cgc_genes)

    output_data = pandas.DataFrame(data=[["" for x in args.annotation] for y in every_genes], index=every_genes, columns=args.annotation, dtype=str)
    for annotation in tqdm.tqdm(args.annotation):
        output_data.loc[input_data[annotation] & set(every_genes), annotation] = "*"
    output_data.index.name = index_name

    print(output_data)
    output_data.to_latex(args.output + ".tex", column_format="l" + "c" * len(args.annotation))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig = matplotlib.pyplot.figure(figsize=(2 ** len(input_data) + 40, 24))

    try:
        upsetplot.plot(upsetplot.from_contents(input_data), fig=fig, show_counts="%d", show_percentages=True, element_size=None)
    except IndexError:
        matplotlib.pyplot.text(0.5, 0.5, "Nothing to show...", fontsize=step00.matplotlib_parameters["axes.titlesize"], color="k", horizontalalignment="center", verticalalignment="center")
        matplotlib.pyplot.xticks([])
        matplotlib.pyplot.yticks([])

    fig.savefig(args.output + ".pdf", bbox_inches="tight")
    matplotlib.pyplot.close(fig)
