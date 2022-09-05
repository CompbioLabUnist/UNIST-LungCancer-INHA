"""
plot_gene_loci_MutationSharedProportion.py: plot gene loci with Mutation Shared Proportion
"""
import argparse
import itertools
import multiprocessing
import tarfile
import typing
from adjustText import adjust_text
import matplotlib
import matplotlib.patches
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import tqdm
import step00

cgc_data = pandas.DataFrame()
mutect_data = pandas.DataFrame()
filtered_mutect_data = pandas.DataFrame()
patients: typing.Set[str] = set()
stage_list: typing.List[str] = list()


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


def run(gene: str, MSP: str) -> str:
    if "SYN" in MSP:
        data = mutect_data.loc[(mutect_data["Hugo_Symbol"] == gene)]
    else:
        data = filtered_mutect_data.loc[(filtered_mutect_data["Hugo_Symbol"] == gene)]

    fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))
    texts = list()

    for i, stage in enumerate(stage_list):
        i += 1

        ax.axhline(i, linestyle="--", color="silver")
        ax.axhline(-1 * i, linestyle="--", color="silver")

        matplotlib.pyplot.text(cgc_data.loc[gene, "start"], i, stage, color="darkred", fontsize="xx-small", horizontalalignment="left", verticalalignment="bottom")
        matplotlib.pyplot.text(cgc_data.loc[gene, "end"], i, stage, color="darkred", fontsize="xx-small", horizontalalignment="right", verticalalignment="bottom")
        matplotlib.pyplot.text(cgc_data.loc[gene, "start"], -1 * i, stage, color="darkblue", fontsize="xx-small", horizontalalignment="left", verticalalignment="top")
        matplotlib.pyplot.text(cgc_data.loc[gene, "end"], -1 * i, stage, color="darkblue", fontsize="xx-small", horizontalalignment="right", verticalalignment="top")

    for patient in patients:
        selected_data = data.loc[(data[MSP] == "Lower") & (data["Patient"] == patient), ["Stage"] + step00.sharing_strategy]
        mutation_set = set(selected_data.loc[(selected_data["Stage"] == stage_list[0]), step00.sharing_strategy].itertuples(index=False, name=None)) & set(selected_data.loc[(selected_data["Stage"] == stage_list[1]), step00.sharing_strategy].itertuples(index=False, name=None))
        for mutation in mutation_set:
            matplotlib.pyplot.axvline(numpy.mean(mutation[1:3]), ymin=0, ymax=0.5, linestyle=":", color="blue", alpha=0.3, zorder=0)

        selected_data = data.loc[(data[MSP] == "Higher") & (data["Patient"] == patient)]
        mutation_set = set(selected_data.loc[(selected_data["Stage"] == stage_list[0]), step00.sharing_strategy].itertuples(index=False, name=None)) & set(selected_data.loc[(selected_data["Stage"] == stage_list[1]), step00.sharing_strategy].itertuples(index=False, name=None))
        for mutation in mutation_set:
            matplotlib.pyplot.axvline(numpy.mean(mutation[1:3]), ymin=0.5, ymax=1.0, linestyle=":", color="red", alpha=0.3, zorder=0)

    for i, stage in enumerate(stage_list):
        i += 1

        drawing_data = data.loc[(data[MSP] == "Lower") & (data["Stage"] == stage)]
        for index, row in drawing_data.iterrows():
            texts.append(matplotlib.pyplot.text(numpy.mean(row[["Start_Position", "End_Position"]]), -1 * i, row["Patient"], color=step00.stage_color_code[row["Stage"]], fontsize="xx-small", horizontalalignment="center", verticalalignment="center"))

        drawing_data = data.loc[(data[MSP] == "Higher") & (data["Stage"] == stage)]
        for index, row in drawing_data.iterrows():
            texts.append(matplotlib.pyplot.text(numpy.mean(row[["Start_Position", "End_Position"]]), i, row["Patient"], color=step00.stage_color_code[row["Stage"]], fontsize="xx-small", horizontalalignment="center", verticalalignment="center"))

    ax.add_patch(matplotlib.patches.Rectangle((cgc_data.loc[gene, "start"], -0.25), (cgc_data.loc[gene, "end"] - cgc_data.loc[gene, "start"]), 0.5, color="silver", zorder=2))
    matplotlib.pyplot.text(numpy.mean(cgc_data.loc[gene, ["start", "end"]]), 0, f"{gene} (chr{cgc_data.loc[gene, 'Genome Location']})", color="k", horizontalalignment="center", verticalalignment="center", zorder=2)

    matplotlib.pyplot.xlim((cgc_data.loc[gene, "start"], cgc_data.loc[gene, "end"]))
    matplotlib.pyplot.ylim((-1 * len(stage_list) - 1, len(stage_list) + 1))
    matplotlib.pyplot.xticks([])
    matplotlib.pyplot.yticks([])
    matplotlib.pyplot.title(MSP)
    matplotlib.pyplot.tight_layout()

    adjust_text(texts, autoalign="y", only_move={"points": "y", "text": "y", "objects": "y"}, arrowprops={"arrowstyle": "-", "color": "k", "linewidth": 0.5}, ax=ax, lim=step00.small)

    figure = f"{gene}-{MSP}.pdf"
    fig.savefig(figure)
    matplotlib.pyplot.close(fig)
    return figure


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("cgc", help="CGC gene CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_strategy = parser.add_mutually_exclusive_group(required=True)
    group_strategy.add_argument("--median", help="Median division", action="store_true", default=False)
    group_strategy.add_argument("--mean", help="Mean division", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.cgc.endswith(".csv"):
        raise ValueError("CGC must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(patients)

    cgc_data = pandas.read_csv(args.cgc, index_col="Gene Symbol")
    cgc_data = cgc_data.loc[~(cgc_data["Genome Location"].str.contains(":-"))]
    cgc_data["chromosome"] = list(map(lambda x: "chr" + (x.replace("-", ":").split(":")[0]), cgc_data["Genome Location"]))
    cgc_data["start"] = list(map(lambda x: int(x.replace("-", ":").split(":")[1]), cgc_data["Genome Location"]))
    cgc_data["end"] = list(map(lambda x: int(x.replace("-", ":").split(":")[2]), cgc_data["Genome Location"]))
    print(cgc_data)

    gene_list = list(cgc_data.index)

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))
    args.input.sort(key=step00.sorting)
    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
        mutect_data["Tumor_Sample_Barcode"] = pool.map(step00.get_id, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Patient"] = pool.map(step00.get_patient, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Stage"] = pool.map(step00.get_long_sample_type, mutect_data["Tumor_Sample_Barcode"])
    print(mutect_data)

    mutect_data = mutect_data[(mutect_data["Hugo_Symbol"].isin(gene_list))]
    print(mutect_data)

    for MSP in tqdm.tqdm(step00.sharing_columns):
        if args.median:
            cutting = numpy.median(clinical_data[MSP])
        elif args.mean:
            cutting = numpy.mean(clinical_data[MSP])
        else:
            raise Exception("Something went wrong!!")

        lower_patients = set(clinical_data.loc[(clinical_data[MSP] <= cutting)].index)
        mutect_data[MSP] = list(map(lambda x: "Lower" if (step00.get_patient(x) in lower_patients) else "Higher", mutect_data["Tumor_Sample_Barcode"]))
    print(mutect_data)

    filtered_mutect_data = mutect_data[(mutect_data[step00.nonsynonymous_column].isin(step00.nonsynonymous_mutations))]
    print(filtered_mutect_data)

    patients &= set(mutect_data["Patient"])
    print(patients)

    stage_list = list(filter(lambda x: x in set(mutect_data["Stage"]), reversed(step00.long_sample_type_list)))
    print(stage_list)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    with multiprocessing.Pool(args.cpus) as pool:
        figures = pool.starmap(run, itertools.product(gene_list, step00.sharing_columns))

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
