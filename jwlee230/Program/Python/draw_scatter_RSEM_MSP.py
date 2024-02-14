"""
draw_scatter_RSEM_MSP.py: draw scatter plots with DEG & MSP
"""
import argparse
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import tqdm
import step00

input_data = pandas.DataFrame()
clinical_data = pandas.DataFrame()


def get_middle(values):
    return (max(values) + min(values)) / 2


def scatter(MSP, gene):
    sample_set = set(clinical_data[f"{MSP}-sample"]) | set(map(step00.get_paired_primary, clinical_data[f"{MSP}-sample"]))
    drawing_data = expression_data.loc[list(filter(lambda x: x in sample_set, list(expression_data.index))), :]

    precancer_r = input_data.loc[gene, f"Precancer-{MSP}-r"]
    precancer_slope = input_data.loc[gene, f"Precancer-{MSP}-slope"] if (precancer_r > 0) else (-1 * input_data.loc[gene, f"Precancer-{MSP}-slope"])

    primary_r = input_data.loc[gene, f"Primary-{MSP}-r"]
    primary_slope = input_data.loc[gene, f"Primary-{MSP}-slope"] if (primary_r > 0) else (-1 * input_data.loc[gene, f"Primary-{MSP}-slope"])

    g = seaborn.lmplot(data=drawing_data, x=MSP, y=gene, hue="Stage", hue_order=["Precancer", "Primary"], palette={"Precancer": "tab:pink", "Primary": "gray"}, scatter=True, fit_reg=True, height=18, aspect=1)
    g.fig.text(0.5, 0.5, f"Precancer: r={precancer_r:.3f}, slope={precancer_slope:.1e}\nPrimary: r={primary_r:.3f}, slope={primary_slope:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
    g.fig.tight_layout()

    fig_name = f"Scatter-{MSP}-{gene}.pdf"
    g.savefig(fig_name)
    matplotlib.pyplot.close(g.fig)
    return fig_name


def joint(MSP, gene):
    sample_set = set(clinical_data[f"{MSP}-sample"]) | set(map(step00.get_paired_primary, clinical_data[f"{MSP}-sample"]))
    drawing_data = expression_data.loc[list(filter(lambda x: x in sample_set, list(expression_data.index))), :]

    precancer_r = input_data.loc[gene, f"Precancer-{MSP}-r"]
    precancer_slope = input_data.loc[gene, f"Precancer-{MSP}-slope"] if (precancer_r > 0) else (-1 * input_data.loc[gene, f"Precancer-{MSP}-slope"])

    primary_r = input_data.loc[gene, f"Primary-{MSP}-r"]
    primary_slope = input_data.loc[gene, f"Primary-{MSP}-slope"] if (primary_r > 0) else (-1 * input_data.loc[gene, f"Primary-{MSP}-slope"])

    g = seaborn.jointplot(data=drawing_data, x=MSP, y=gene, hue="Stage", hue_order=["Precancer", "Primary"], palette={"Precancer": "tab:pink", "Primary": "gray"}, height=18, ratio=5)
    g.fig.text(0.5, 0.5, f"Precancer: r={precancer_r:.3f}, slope={precancer_slope:.1e}\nPrimary: r={primary_r:.3f}, slope={primary_slope:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
    g.fig.tight_layout()

    fig_name = f"Joint-{MSP}-{gene}.pdf"
    g.savefig(fig_name)
    matplotlib.pyplot.close(g.fig)
    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input DEG-MSP TSV file", type=str)
    parser.add_argument("expression", help="Expression TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--r", help="r-value threshold", type=float, default=0.3)
    parser.add_argument("--slope", help="Slope threshold", type=float, default=5)
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

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    patients = set(clinical_data.index)
    print(clinical_data)

    expression_data = pandas.read_csv(args.expression, sep="\t", index_col=0).T
    expression_data = expression_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(expression_data.index))), :]
    expression_data["Stage"] = list(map(lambda x: "Primary" if (step00.get_long_sample_type(x) == "Primary") else "Precancer", list(expression_data.index)))
    for column in tqdm.tqdm(step00.sharing_columns):
        expression_data[column] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), column], list(expression_data.index)))
    print(expression_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    with multiprocessing.Pool(processes=args.cpus) as pool:
        for MSP in tqdm.tqdm(step00.sharing_columns[:1]):
            primary_POS_gene_set = set(input_data.loc[(input_data[f"Primary-{MSP}-slope"] > args.slope) & (input_data[f"Primary-{MSP}-r"] > args.r)])
            primary_NEG_gene_set = set(input_data.loc[(input_data[f"Primary-{MSP}-slope"] > args.slope) & (input_data[f"Primary-{MSP}-r"] < (-1 * args.r))])

            genes = list()
            genes += sorted(set(input_data.loc[(input_data[f"Precancer-{MSP}-slope"] > args.slope) & (input_data[f"Precancer-{MSP}-r"] > args.r)].index) - primary_POS_gene_set)
            genes += sorted(set(input_data.loc[(input_data[f"Precancer-{MSP}-slope"] > args.slope) & (input_data[f"Precancer-{MSP}-r"] < (-1 * args.r))].index) - primary_NEG_gene_set)

            figures += list(pool.starmap(scatter, [(MSP, gene) for gene in genes]))
            figures += list(pool.starmap(joint, [(MSP, gene) for gene in genes]))

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
