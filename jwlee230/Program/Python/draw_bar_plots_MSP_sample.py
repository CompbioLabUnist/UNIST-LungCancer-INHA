"""
draw_bar_plots_MSP_sample.py: draw bar plots upon DEG with MSP with precancer vs. primary
"""
import argparse
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import statannotations.Annotator
import tqdm
import step00

expression_data = pandas.DataFrame()
metadata = pandas.DataFrame()


def run(MSP: str, gene: str) -> str:
    lower_bound, higher_bound = numpy.quantile(clinical_data[MSP], args.percentage), numpy.quantile(clinical_data[MSP], 1.0 - args.percentage)

    expression_samples = set(expression_data.index)

    lower_precancer_list = list(filter(lambda x: x in expression_samples, list(clinical_data.loc[(clinical_data[MSP] <= lower_bound), f"{MSP}-sample"])))
    higher_precancer_list = list(filter(lambda x: x in expression_samples, list(clinical_data.loc[(clinical_data[MSP] >= higher_bound), f"{MSP}-sample"])))

    lower_normal_list = list(filter(lambda x: x in expression_samples, list(map(step00.get_paired_normal, lower_precancer_list))))
    higher_normal_list = list(filter(lambda x: x in expression_samples, list(map(step00.get_paired_normal, higher_precancer_list))))

    lower_primary_list = list(filter(lambda x: x in expression_samples, list(map(step00.get_paired_primary, lower_precancer_list))))
    higher_primary_list = list(filter(lambda x: x in expression_samples, list(map(step00.get_paired_primary, higher_precancer_list))))

    raw_output_data = list()
    raw_output_data += [(sample, "Lower", "Normal", expression_data.loc[sample, gene]) for sample in lower_normal_list]
    raw_output_data += [(sample, "Lower", "Precancer", expression_data.loc[sample, gene]) for sample in lower_precancer_list]
    raw_output_data += [(sample, "Lower", "Primary", expression_data.loc[sample, gene]) for sample in lower_primary_list]
    raw_output_data += [(sample, "Higher", "Normal", expression_data.loc[sample, gene]) for sample in higher_normal_list]
    raw_output_data += [(sample, "Higher", "Precancer", expression_data.loc[sample, gene]) for sample in higher_precancer_list]
    raw_output_data += [(sample, "Higher", "Primary", expression_data.loc[sample, gene]) for sample in higher_primary_list]

    output_data = pandas.DataFrame(raw_output_data, columns=["Sample", "Lower/Higher", "PRE/PRI", "Expression"])

    p1 = scipy.stats.mannwhitneyu(output_data.loc[(output_data["Lower/Higher"] == "Lower") & (output_data["PRE/PRI"] == "Precancer"), "Expression"], output_data.loc[(output_data["Lower/Higher"] == "Lower") & (output_data["PRE/PRI"] == "Primary"), "Expression"])[1]
    p2 = scipy.stats.mannwhitneyu(output_data.loc[(output_data["Lower/Higher"] == "Lower") & (output_data["PRE/PRI"] == "Precancer"), "Expression"], output_data.loc[(output_data["Lower/Higher"] == "Higher") & (output_data["PRE/PRI"] == "Precancer"), "Expression"])[1]
    p3 = scipy.stats.mannwhitneyu(output_data.loc[(output_data["Lower/Higher"] == "Higher") & (output_data["PRE/PRI"] == "Precancer"), "Expression"], output_data.loc[(output_data["Lower/Higher"] == "Higher") & (output_data["PRE/PRI"] == "Primary"), "Expression"])[1]
    p4 = scipy.stats.mannwhitneyu(output_data.loc[(output_data["Lower/Higher"] == "Lower") & (output_data["PRE/PRI"] == "Primary"), "Expression"], output_data.loc[(output_data["Lower/Higher"] == "Higher") & (output_data["PRE/PRI"] == "Primary"), "Expression"])[1]

    if (p1 >= 0.05) or (p2 >= 0.05) or (p3 < 0.05) or (p4 < 0.05):
        return ""

    print("Lower-NML", len(lower_normal_list))
    print("Higher-NML", len(higher_normal_list))
    print("Lower-PRE:", len(lower_precancer_list))
    print("Higher-PRE:", len(higher_precancer_list))
    print("Lower-PRI:", len(lower_primary_list))
    print("Higher-PRI:", len(higher_primary_list))

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    seaborn.violinplot(data=output_data, x="Lower/Higher", y="Expression", hue="PRE/PRI", order=["Lower", "Higher"], hue_order=["Normal", "Precancer", "Primary"], palette={"Normal": "cyan", "Precancer": "tab:pink", "Primary": "gray"}, inner="box", linewidth=5, cut=1, ax=ax)
    statannotations.Annotator.Annotator(ax, [(("Lower", "Normal"), ("Lower", "Precancer")), (("Higher", "Normal"), ("Higher", "Precancer")), (("Lower", "Normal"), ("Higher", "Normal"))] + list(filter(lambda x: (x[0][0] == x[1][0]) or (x[0][1] == x[1][1]), itertools.combinations(itertools.product(["Lower", "Higher"], ["Precancer", "Primary"]), r=2))), data=output_data, x="Lower/Higher", y="Expression", hue="PRE/PRI", order=["Lower", "Higher"], hue_order=["Normal", "Precancer", "Primary"]).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()

    matplotlib.pyplot.ylabel("Gene expression")
    matplotlib.pyplot.title(gene)
    matplotlib.pyplot.legend(loc="upper left")
    matplotlib.pyplot.tight_layout()

    fig_name = f"{MSP}-{gene}.pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input DEG-MSP TSV file", type=str)
    parser.add_argument("expression", help="Expression TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--r", help="r-value threshold", type=float, default=0.4)
    parser.add_argument("--slope", help="Slope threshold", type=float, default=7.5)
    parser.add_argument("--percentage", help="Percentage of patients to include", type=float, default=0.25)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

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
    elif not (0.0 < args.percentage < 0.5):
        raise ValueError("Percentage must be (0.0, 0.5)!!")
    elif args.cpus < 1:
        raise ValueError("Number of CPUs must be positive!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(len(patients), sorted(patients))

    expression_data = pandas.read_csv(args.expression, sep="\t", index_col=0).T
    expression_data = expression_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(expression_data.index))), :]
    expression_data["Stage"] = list(map(step00.get_long_sample_type, list(expression_data.index)))
    for column in tqdm.tqdm(step00.sharing_columns):
        expression_data[column] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), column], list(expression_data.index)))
    print(expression_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    with multiprocessing.Pool(args.cpus) as pool:
        for MSP in tqdm.tqdm(step00.sharing_columns[:1]):
            primary_POS_gene_set = set(input_data.loc[(input_data[f"Primary-{MSP}-slope"] > args.slope) & (input_data[f"Primary-{MSP}-r"] > args.r)])
            primary_NEG_gene_set = set(input_data.loc[(input_data[f"Primary-{MSP}-slope"] > args.slope) & (input_data[f"Primary-{MSP}-r"] < (-1 * args.r))])

            genes = list()
            genes += sorted(set(input_data.loc[(input_data[f"Precancer-{MSP}-slope"] > args.slope) & (input_data[f"Precancer-{MSP}-r"] > args.r)].index) - primary_POS_gene_set)
            genes += sorted(set(input_data.loc[(input_data[f"Precancer-{MSP}-slope"] > args.slope) & (input_data[f"Precancer-{MSP}-r"] < (-1 * args.r))].index) - primary_NEG_gene_set)

            figures += list(filter(None, pool.starmap(run, [(MSP, gene) for gene in genes])))

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(figures):
            tar.add(f, arcname=f)
