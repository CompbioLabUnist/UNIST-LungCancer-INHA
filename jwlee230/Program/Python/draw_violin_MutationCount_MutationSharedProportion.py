"""
draw_violin_MutationCount_MutationSharedProportion.py: draw violin plot about mutation count with Mutation Shared Proportion
"""
import argparse
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import statannotations.Annotator
import tqdm
import step00


def get_count(filename):
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False).shape[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 MAF file(s)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical information w/ Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_threshold = parser.add_mutually_exclusive_group(required=True)
    group_threshold.add_argument("--median", help="Use Median threshold", action="store_true", default=False)
    group_threshold.add_argument("--mean", help="Use Mean threshold", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("Input must end with .MAF!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    print(clinical_data)

    args.input = list(filter(lambda x: step00.get_patient(x) in set(clinical_data.index), args.input))

    input_data = pandas.DataFrame(index=list(map(lambda x: step00.get_id(x), args.input)))
    with multiprocessing.Pool(args.cpus) as pool:
        input_data["Count"] = pool.map(get_count, args.input)
        input_data["Patient"] = pool.map(step00.get_patient, list(input_data.index))
        input_data["Stage"] = pool.map(step00.get_long_sample_type, list(input_data.index))
    print(input_data)

    for MSP in tqdm.tqdm(step00.sharing_columns):
        if args.median:
            threshold = numpy.median(clinical_data[MSP])
        elif args.mean:
            threshold = numpy.mean(clinical_data[MSP])
        else:
            raise Exception("Something went wrong!!")
        input_data[MSP] = list(map(lambda x: "Lower" if (clinical_data.loc[step00.get_patient(x), MSP] < threshold) else "Higher", list(input_data.index)))
    print(input_data)

    figures = list()
    for MSP in tqdm.tqdm(step00.sharing_columns):
        stage_list = list(filter(lambda x: (x in set(input_data.loc[(input_data[MSP] == "Lower"), "Stage"])) and (x in set(input_data.loc[(input_data[MSP] == "Higher"), "Stage"])), step00.long_sample_type_list))

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        seaborn.violinplot(data=input_data, x="Stage", y="Count", order=stage_list, hue=MSP, hue_order=["Lower", "Higher"], cut=1, linewidth=5, ax=ax)
        statannotations.Annotator.Annotator(ax, list(map(lambda x: ((x, "Lower"), (x, "Higher")), stage_list)), data=input_data, x="Stage", y="Count", order=stage_list, hue=MSP, hue_order=["Lower", "Higher"]).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

        matplotlib.pyplot.ylabel("Short Nucleotide Variation Counts")
        matplotlib.pyplot.tight_layout()

        figures.append(f"{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
