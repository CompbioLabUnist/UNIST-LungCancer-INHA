"""
draw_violin_plots_geneset.py: draw violin plots upon DEG with given gene set
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
    stage_order = list(filter(lambda x: x in set(input_data["Stage"]), step00.long_sample_type_list))

    compare_list = list()
    for s1, s2 in itertools.combinations(stage_order, r=2):
        stat, p = scipy.stats.mannwhitneyu(input_data.loc[(input_data["Stage"] == s1), gene], input_data.loc[(input_data["Stage"] == s2), gene])
        if p < 0.05:
            compare_list.append((s1, s2))

    try:
        stat, p = scipy.stats.kruskal(*[input_data.loc[(input_data["Stage"] == stage), gene] for stage in stage_order])
    except ValueError:
        _, p = 0., 1.0

    if (p > 0.05):
        return ""

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    seaborn.violinplot(data=input_data, x="Stage", order=stage_order, palette=step00.stage_color_code, y=gene, cut=1, linewidth=5, ax=ax)
    if compare_list:
        statannotations.Annotator.Annotator(ax, compare_list, data=input_data, x="Stage", y=gene, order=stage_order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()

    matplotlib.pyplot.ylabel(f"{gene} expression")
    matplotlib.pyplot.title(f"{gene}: K.W. p={p:.3f}")
    matplotlib.pyplot.tight_layout()

    fig_name = f"{gene}.pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TPM TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data data CSV file", type=str)
    parser.add_argument("geneset", help="Gene set TXT file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .tsv!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.geneset.endswith(".txt"):
        raise ValueError("Gene set must end with .TXT!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
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

    geneset_data = pandas.read_csv(args.geneset, comment="#")
    gene_set = set(geneset_data[list(geneset_data.columns)[0]])

    input_data = pandas.read_csv(args.input, sep="\t", index_col="gene_name").T
    gene_set &= set(input_data.columns)
    input_data = input_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index))), sorted(gene_set)]
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    input_data["All markers"] = list(map(lambda x: sum(input_data.loc[x, sorted(gene_set)]), list(input_data.index)))
    print(input_data)

    with multiprocessing.Pool(args.cpus) as pool:
        tar_files = list(filter(None, pool.map(run, sorted(gene_set))))

    stage_order = list(filter(lambda x: x in set(input_data["Stage"]), step00.long_sample_type_list))
    gene = "All markers"

    compare_list = list()
    for s1, s2 in itertools.combinations(stage_order, r=2):
        stat, p = scipy.stats.mannwhitneyu(input_data.loc[(input_data["Stage"] == s1), gene], input_data.loc[(input_data["Stage"] == s2), gene])
        if p < 0.05:
            compare_list.append((s1, s2))

    try:
        stat, p = scipy.stats.kruskal(*[input_data.loc[(input_data["Stage"] == stage), gene] for stage in stage_order])
    except ValueError:
        _, p = 0., 1.0

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    seaborn.violinplot(data=input_data, x="Stage", order=stage_order, palette=step00.stage_color_code, y=gene, cut=1, linewidth=5, ax=ax)
    if compare_list:
        statannotations.Annotator.Annotator(ax, compare_list, data=input_data, x="Stage", y=gene, order=stage_order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()

    matplotlib.pyplot.ylabel(f"{gene} expression")
    matplotlib.pyplot.title(f"{gene}: K.W. p={p:.3f}")
    matplotlib.pyplot.tight_layout()

    fig_name = f"{gene}.pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)
    tar_files.append(fig_name)

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(tar_files):
            tar.add(f, arcname=f)
