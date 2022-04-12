"""
draw_violin_plots.py: draw violin plots upon DEG
"""
import argparse
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
ylabel = ""


def run(gene: str) -> str:
    stage_order = list(filter(lambda x: x in set(input_data["Stage"]), step00.long_sample_type_list))

    stat, p = scipy.stats.kruskal(*[input_data.loc[(input_data["Stage"] == stage), cell] for stage in stage_order])

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.violinplot(data=input_data, x="Stage", y=gene, order=stage_order, ax=ax)
    statannotations.Annotator.Annotator(ax, list(zip(stage_order, stage_order[1:])), data=input_data, x="Stage", y=gene, order=stage_order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.ylabel(ylabel)
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

    ylabel = args.input.split("/")[-1].split(".")[1]

    input_data = pandas.read_csv(args.input, sep="\t", index_col="gene_name").T
    gene_set &= set(input_data.columns)
    input_data = input_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index))), sorted(gene_set)]
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    print(input_data)

    with multiprocessing.Pool(args.cpus) as pool:
        tar_files = pool.map(run, sorted(gene_set))

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(tar_files):
            tar.add(f, arcname=f)
