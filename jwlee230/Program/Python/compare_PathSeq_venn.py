"""
compare_PathSeq_venn.py: compare Pathseq WES vs. WTS with venn diagram
"""
import argparse
import matplotlib
import matplotlib.pyplot
import pandas
import tqdm
import upsetplot
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PathSeq results TSV file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--annotation", help="Annotation for figure", type=str, nargs="+", required=True)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif len(args.input) != len(args.annotation):
        raise ValueError("INPUT must have same ANNOTATION!!")

    output_data = dict()
    for input_file, input_name in tqdm.tqdm(zip(args.input, args.annotation)):
        input_data = pandas.read_csv(input_file, sep="\t", index_col=0)
        output_data[input_name] = set(list(input_data.columns)[:-1])

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig = matplotlib.pyplot.figure(figsize=(2 ** len(args.input) + 40, 24))

    upsetplot.plot(upsetplot.from_contents(output_data), fig=fig, show_counts="%d", show_percentages=True, element_size=None)

    fig.savefig(args.output, bbox_inches="tight")
    matplotlib.pyplot.close(fig)
