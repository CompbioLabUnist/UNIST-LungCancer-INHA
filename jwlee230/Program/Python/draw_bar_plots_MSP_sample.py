"""
draw_bar_plots_MSP_sample.py: draw bar plots upon DEG with MSP with precancer vs. primary
"""
import argparse
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import tqdm
import step00

expression_data = pandas.DataFrame()
order = ["Lower", "Higher"]


def run(sample: str, gene: str) -> str:
    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    matplotlib.pyplot.bar(0, expression_data.loc[sample, gene], width=0.8, color=step00.stage_color_code[step00.get_long_sample_type(sample)], label=step00.get_long_sample_type(sample))
    matplotlib.pyplot.bar(1, expression_data.loc[step00.get_paired_primary(sample), gene], width=0.8, color=step00.stage_color_code["Primary"], label="Primary")

    matplotlib.pyplot.xticks([0, 1], [step00.get_long_sample_type(sample), "Primary"])
    matplotlib.pyplot.xlabel(f"{step00.get_patient(sample)}")
    matplotlib.pyplot.ylabel(f"{gene} expression")
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.legend()
    matplotlib.pyplot.tight_layout()

    fig_name = f"{sample}-{gene}.pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input DEG-MSP TSV file", type=str)
    parser.add_argument("expression", help="Expression TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--r", help="r-value threshold", type=float, default=0.7)
    parser.add_argument("--slope", help="Slope threshold", type=float, default=100)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

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
    elif args.cpus < 1:
        raise ValueError("Number of CPUs must be positive!!")

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

    genes = set()
    for MSP in tqdm.tqdm(step00.sharing_columns):
        genes |= set(input_data.loc[(input_data[f"CIS+AIS-{MSP}-slope"] > args.slope) & ((input_data[f"CIS+AIS-{MSP}-r"] > args.r) | (input_data[f"CIS+AIS-{MSP}-r"] < (-1 * args.r)))].index)
    print(len(genes))

    samples = list(expression_data.index)
    precancer_samples = list(filter(lambda x: (step00.get_long_sample_type(x) != "Primary") and (step00.get_paired_primary(x) in samples), samples))
    print(len(precancer_samples))

    figures = list()
    with multiprocessing.Pool(args.cpus) as pool:
        figures = list(pool.starmap(run, itertools.product(precancer_samples, sorted(genes))))

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(figures):
            tar.add(f, arcname=f)
