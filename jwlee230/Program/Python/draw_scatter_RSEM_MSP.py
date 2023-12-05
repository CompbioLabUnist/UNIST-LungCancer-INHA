"""
draw_scatter_RSEM_MSP.py: draw scatter plots with DEG & MSP
"""
import argparse
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import tqdm
import step00

input_data = pandas.DataFrame()
expression_data = pandas.DataFrame()


def get_middle(values):
    return (max(values) + min(values)) / 2


def scatter(stage, MSP, gene):
    if stage == "All":
        tmp_data = expression_data
        color = "blue"
    elif stage == "Precancer":
        tmp_data = expression_data[~(expression_data["Stage"].isin({"Normal", "Primary"}))]
        color = "tab:pink"
    else:
        tmp_data = expression_data[(expression_data["Stage"] == stage)]
        color = step00.stage_color_code[stage]

    r = input_data.loc[gene, f"{stage}-{MSP}-r"]
    slope = input_data.loc[gene, f"{stage}-{MSP}-slope"] if (r > 0) else (-1 * input_data.loc[gene, f"{stage}-{MSP}-slope"])

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    seaborn.regplot(data=tmp_data, x=MSP, y=gene, color=color, ax=ax)
    matplotlib.pyplot.text(x=get_middle(tmp_data[MSP]), y=get_middle(tmp_data[gene]), s=f"r={r:.3f}, slope={slope:.1e}", horizontalalignment="center", verticalalignment="center", fontsize="small", bbox={"alpha": 0.3, "color": "white"})

    matplotlib.pyplot.tight_layout()

    fig_name = f"Scatter-{stage}-{MSP}-{gene}.pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)
    return fig_name


def joint(stage, MSP, gene):
    if stage == "All":
        tmp_data = expression_data
        color = "blue"
    elif stage == "Precancer":
        tmp_data = expression_data[~(expression_data["Stage"].isin({"Normal", "Primary"}))]
        color = "tab:pink"
    else:
        tmp_data = expression_data[(expression_data["Stage"] == stage)]
        color = step00.stage_color_code[stage]

    r = input_data.loc[gene, f"{stage}-{MSP}-r"]
    slope = input_data.loc[gene, f"{stage}-{MSP}-slope"] if (r > 0) else (-1 * input_data.loc[gene, f"{stage}-{MSP}-slope"])

    g = seaborn.jointplot(data=tmp_data, x=MSP, y=gene, color=color, kind="reg", height=24, ratio=5)
    g.fig.text(0.5, 0.5, f"r={r:.3f}, slope={slope:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})

    fig_name = f"Joint-{stage}-{MSP}-{gene}.pdf"
    g.savefig(fig_name)
    matplotlib.pyplot.close(g.fig)
    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input DEG-MSP TSV file", type=str)
    parser.add_argument("expression", help="Expression TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--r", help="r-value threshold", type=float, default=0.5)
    parser.add_argument("--slope", help="Slope threshold", type=float, default=10)
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

    stages = step00.long_sample_type_list + ["Precancer", "All"]

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    with multiprocessing.Pool(processes=args.cpus) as pool:
        for stage, MSP in tqdm.tqdm(list(itertools.product(stages, step00.sharing_columns[:1]))):
            if f"{stage}-{MSP}-log10(abs(slope))" not in set(input_data.columns):
                continue

            genes = list(input_data.loc[(input_data[f"{stage}-{MSP}-slope"] > args.slope) & ((input_data[f"{stage}-{MSP}-r"] > args.r) | (input_data[f"{stage}-{MSP}-r"] < (-1 * args.r)))].index)

            figures += list(pool.starmap(scatter, [(stage, MSP, gene) for gene in genes]))
            figures += list(pool.starmap(joint, [(stage, MSP, gene) for gene in genes]))

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
