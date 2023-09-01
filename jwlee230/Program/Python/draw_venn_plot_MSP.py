"""
draw_venn_plot_MSP.py: draw a venn diagram between DEG and MSP
"""
import argparse
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import upsetplot
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input DEG-MSP TSV file", type=str)
    parser.add_argument("expression", help="Expression TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--r", help="r-value threshold", type=float, default=0.3)
    parser.add_argument("--slope", help="Slope threshold", type=float, default=5)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.expression.endswith(".tsv"):
        raise ValueError("Expression must end with .TSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif not (0 < args.r < 1):
        raise ValueError("r-value must be in (0, 1)!!")
    elif args.slope <= 0:
        raise ValueError("Slope must be positive!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)

    clinical_data: pandas.DataFrame = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    patients = set(clinical_data.index)
    print(clinical_data)

    expression_data = pandas.read_csv(args.expression, sep="\t", index_col=0).T
    expression_data = expression_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(expression_data.index))), :]
    expression_data["Stage"] = list(map(step00.get_long_sample_type, list(expression_data.index)))
    for column in tqdm.tqdm(step00.sharing_columns):
        expression_data[column] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), column], list(expression_data.index)))
    print(expression_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    figures = list()
    for MSP in tqdm.tqdm(step00.sharing_columns):
        venn_data = dict()
        for stage in ["CIS+AIS", "Primary"]:
            if f"{stage}-{MSP}-log10(abs(slope))" not in set(input_data.columns):
                continue

            venn_data[f"{stage}-pos"] = set(input_data.loc[(input_data[f"{stage}-{MSP}-slope"] > args.slope) & (input_data[f"{stage}-{MSP}-r"] > args.r)].index)
            venn_data[f"{stage}-neg"] = set(input_data.loc[(input_data[f"{stage}-{MSP}-slope"] > args.slope) & (input_data[f"{stage}-{MSP}-r"] < (-1 * args.r))].index)

        fig = matplotlib.pyplot.figure(figsize=(10 * len(venn_data) + 10, 24))

        try:
            upsetplot.plot(upsetplot.from_contents(venn_data), fig=fig, show_counts="%d", show_percentages=True, element_size=None)
        except IndexError:
            pass

        figures.append(f"{MSP}.pdf")
        fig.savefig(figures[-1], bbox_inches="tight")
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
