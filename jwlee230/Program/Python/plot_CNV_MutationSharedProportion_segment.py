"""
plot_CNV_MutationSharedProportion_segment.py: violin plot Copy Number Variation with Mutation Shared Proportion for segment count
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

input_data = pandas.DataFrame()
watching = ""
threshold = 0.0


def get_chromosome_data(sample: str) -> int:
    tmp_data = input_data.loc[(input_data["Sample"] == sample) & (((input_data[watching] * input_data["weight"]) <= (1 - threshold)) | ((input_data[watching] * input_data["weight"]) >= (1 + args.threshold))), :]
    return tmp_data.shape[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="CNV segment.tsv file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--watching", help="Watching column name", type=str, required=True)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)
    parser.add_argument("--threshold", help="Threshold for gain/loss", type=float, default=0.2)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_strategy = parser.add_mutually_exclusive_group(required=True)
    group_strategy.add_argument("--median", help="Median division", action="store_true", default=False)
    group_strategy.add_argument("--mean", help="Mean division", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.threshold < 1):
        raise ValueError("Threshold must be (0, 1)")

    watching = args.watching
    threshold = args.threshold

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data = input_data.loc[(input_data["Patient"].isin(patients))]
    print(input_data)

    sample_list = sorted(set(input_data["Sample"]), key=step00.sorting_by_type)
    print(sample_list)

    output_data = pandas.DataFrame(data=sample_list, columns=["Sample"])
    with multiprocessing.Pool(args.cpus) as pool:
        output_data["Patient"] = pool.map(step00.get_patient, output_data["Sample"])
        output_data["Stage"] = pool.map(step00.get_long_sample_type, output_data["Sample"])
        output_data["Segment"] = pool.map(get_chromosome_data, output_data["Sample"])
    print(output_data)

    sample_list = sorted(set(output_data["Sample"]), key=step00.sorting_by_type)
    print(sample_list)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()

    for MSP in tqdm.tqdm(step00.sharing_columns):
        if args.median:
            cutting = numpy.median(clinical_data[MSP])
        elif args.mean:
            cutting = numpy.mean(clinical_data[MSP])
        else:
            raise Exception("Something went wrong!!")

        output_data[MSP] = list(map(lambda x: "Lower" if (clinical_data.loc[x, MSP] < cutting) else "Higher", output_data["Patient"]))

        stage_list = list(filter(lambda x: (x in set(output_data.loc[(output_data[MSP] == "Lower"), "Stage"])) and (x in set(output_data.loc[(output_data[MSP] == "Higher"), "Stage"])), step00.long_sample_type_list))
        palette = dict([(stage, step00.stage_color_code[stage]) for stage in stage_list])

        fig, ax = matplotlib.pyplot.subplots(figsize=(20, 18))

        seaborn.violinplot(data=output_data, x=MSP, order=["Lower", "Higher"], y="Segment", hue="Stage", hue_order=stage_list, palette=palette, inner="box", cut=1, ax=ax)

        statannotations.Annotator.Annotator(ax, [(("Lower", stage), ("Higher", stage)) for stage in stage_list], data=output_data, x=MSP, order=["Lower", "Higher"], y="Segment", hue="Stage", hue_order=stage_list).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

        matplotlib.pyplot.title(MSP)
        matplotlib.pyplot.tight_layout()

        figures.append(f"{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
