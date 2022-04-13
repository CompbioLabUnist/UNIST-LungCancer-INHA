"""
draw_MiXCR_diversity.py: draw alpha-diversity from MiXCR results
"""
import argparse
import itertools
import multiprocessing
import tarfile
import typing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import skbio
import statannotations.Annotator
import tqdm
import step00

sample_list: typing.List[str] = list()
input_counts: typing.List[typing.List[int]] = list(list())


def read_result(filename: str) -> pandas.DataFrame:
    data = pandas.read_csv(filename, sep="\t")
    data = data.loc[(data["allVHitsWithScore"].str.contains("TRB"))]
    data["sample"] = step00.get_id(filename)
    data["Stage"] = step00.get_long_sample_type(filename)
    return data


def get_alpha(alpha: str) -> str:
    output_data = pandas.DataFrame(index=sample_list)
    output_data["Stage"] = list(map(step00.get_long_sample_type, sample_list))
    output_data[alpha] = skbio.diversity.alpha_diversity(alpha, input_counts, ids=sample_list)
    output_data = output_data[~(output_data.isin([numpy.nan, numpy.inf, -numpy.inf]).any(1))]

    order = list(filter(lambda x: x in list(output_data["Stage"]), step00.long_sample_type_list))

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    try:
        stat, p = scipy.stats.kruskal(*[output_data.loc[(output_data["Stage"] == stage), alpha] for stage in order])
    except ValueError:
        _, p = 0.0, 1.0

    seaborn.violinplot(data=output_data, x="Stage", y=alpha, order=order, palette=step00.stage_color_code, ax=ax)
    statannotations.Annotator.Annotator(ax, list(itertools.combinations(order, r=2)), data=output_data, x="Stage", y=alpha, order=order, palette=step00.stage_color_code).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.title(f"Kruskal-Wallis p={p:.3f}")
    matplotlib.pyplot.tight_layout()

    figname = f"{alpha}.pdf"
    fig.savefig(figname)
    matplotlib.pyplot.close(fig)

    return figname


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="MiXCR result TXT file", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".txt"), args.input)):
        raise ValueError("Input must end with .TXT!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(pool.map(read_result, args.input), join="outer", ignore_index=True)
    print(input_data)

    sample_list = sorted(set(input_data["sample"]))
    input_counts = [list(input_data.loc[(input_data["sample"] == sample), "cloneCount"].astype(int)) for sample in sample_list]
    max_length = max(list(map(len, input_counts)))
    for i in tqdm.trange(len(sample_list)):
        input_counts[i] = input_counts[i] + [0 for _ in range(max_length - len(input_counts[i]))]
    print(sample_list)

    alphas = ["berger_parker_d", "brillouin_d", "chao1", "dominance", "doubles", "enspie", "fisher_alpha", "gini_index", "goods_coverage", "heip_e", "kempton_taylor_q", "lladser_pe", "margalef", "mcintosh_d", "mcintosh_e", "menhinick", "michaelis_menten_fit", "observed_otus", "pielou_e", "robbins", "shannon", "simpson", "simpson_e", "singles", "strong"]
    print(alphas)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    with multiprocessing.Pool(args.cpus) as pool:
        figures = pool.map(get_alpha, alphas)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
