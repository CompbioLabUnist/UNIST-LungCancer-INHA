"""
draw_violin_plots_clinical.py: draw violin plots upon DEG with clinical data
"""
import argparse
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import scipy.stats
import seaborn
import statannotations.Annotator
import tqdm
import step00

input_data = pandas.DataFrame()


def run(gene: str) -> str:
    compare_list = list()
    for (c1, s1), (c2, s2) in [((a, stage), (b, stage)) for a, b in itertools.combinations(args.compare[1:], r=2) for stage in stage_order] + [((compare, a), (compare, b)) for a, b in itertools.combinations(stage_order, r=2) for compare in args.compare[1:]]:
        stat, p = scipy.stats.mannwhitneyu(input_data.loc[(input_data["Stage"] == s1) & (input_data[args.compare[0]] == c1), gene], input_data.loc[(input_data["Stage"] == s2) & (input_data[args.compare[0]] == c2), gene])
        if p < 0.05:
            compare_list.append(((c1, s1), (c2, s2)))

    try:
        stat, p = scipy.stats.kruskal(*[input_data.loc[(input_data["Stage"] == stage), gene] for stage in stage_order])
    except ValueError:
        _, p = 0.0, 1.0

    if (p > 0.05):
        return ""

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.violinplot(data=input_data, x=args.compare[0], order=args.compare[1:], y=gene, hue="Stage", hue_order=stage_order, palette=step00.stage_color_code, cut=1, linewidth=5, ax=ax)
    if compare_list:
        statannotations.Annotator.Annotator(ax, compare_list, data=input_data, x=args.compare[0], order=args.compare[1:], y=gene, hue="Stage", hue_order=stage_order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.ylabel(f"{gene} expression")
    matplotlib.pyplot.title(f"{gene}: Kruskal-Wallis p={p:.3f}")
    matplotlib.pyplot.tight_layout()

    fig_name = gene + ".pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TPM TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data data CSV file", type=str)
    parser.add_argument("cgc", help="CGC CSV files", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--compare", help="Comparison grouping (type, controls, ...)", type=str, nargs="+", default=["Recurrence", "NO", "YES"])
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .tsv!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.cgc.endswith(".csv"):
        raise ValueError("CGC must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif len(args.compare) < 3:
        raise ValueError("Too few compare values!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    cgc_data = pandas.read_csv(args.cgc, index_col="Gene Symbol")
    gene_set = set(cgc_data.index)
    print(cgc_data)

    input_data = pandas.read_csv(args.input, sep="\t", index_col="gene_name").T

    samples = list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index)))
    gene_set &= set(input_data.columns)

    input_data = input_data.loc[samples, sorted(gene_set)]
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    input_data[args.compare[0]] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), args.compare[0]], list(input_data.index)))
    print(input_data)

    stage_order = list(filter(lambda x: all([(x in set(input_data.loc[(input_data[args.compare[0]] == compare), "Stage"])) for compare in args.compare[1:]]), step00.long_sample_type_list))
    print(stage_order)

    with multiprocessing.Pool(args.cpus) as pool:
        figures = list(filter(None, pool.map(run, sorted(gene_set))))

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(figures):
            tar.add(f, arcname=f)
