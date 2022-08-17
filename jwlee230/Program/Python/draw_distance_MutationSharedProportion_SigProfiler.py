"""
draw_distance_MutationSharedProportion_SigProfiler.py: draw violin plot upon distance of Signature with Mutation Shared Proportion
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
import seaborn
import statannotations.Annotator
import tqdm
import step00

input_data = pandas.DataFrame()
signatures: typing.List[str] = list()


def L1(index):
    return sum(list(map(lambda x: abs(x[0] - x[1]), zip(input_data.loc[index, signatures], input_data.loc[step00.get_paired_primary(index), signatures]))))


def L2(index):
    return numpy.sqrt(sum(list(map(lambda x: (x[0] - x[1]) ** 2, zip(input_data.loc[index, signatures], input_data.loc[step00.get_paired_primary(index), signatures])))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("signature", help="Signature TXT (not necessarily TSV) file", type=str)
    parser.add_argument("clinical", help="Clinical data w/ Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_threshold = parser.add_mutually_exclusive_group(required=True)
    group_threshold.add_argument("--median", help="Use Median threshold", action="store_true", default=False)
    group_threshold.add_argument("--mean", help="Use Mean threshold", action="store_true", default=False)

    group_relative = parser.add_mutually_exclusive_group(required=True)
    group_relative.add_argument("--absolute", help="Use Absolute(Count) value", action="store_true", default=False)
    group_relative.add_argument("--relative", help="Use Relative(Proportion) value", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical data must end with .TSV!!")
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

    input_data = pandas.read_csv(args.signature, sep="\t", index_col=0)
    if args.relative:
        for index in tqdm.tqdm(list(input_data.index)):
            input_data.loc[index, :] = input_data.loc[index, :] / sum(input_data.loc[index, :])
    print(input_data)

    signatures = list(input_data.columns)
    input_data["Patient"] = list(map(step00.get_patient, list(input_data.index)))
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    print(input_data)

    input_data = input_data.loc[(input_data["Patient"].isin(set(clinical_data.index)))]
    print(input_data)

    for MSP in tqdm.tqdm(step00.sharing_columns):
        if args.median:
            threshold = numpy.median(clinical_data[MSP])
        elif args.mean:
            threshold = numpy.mean(clinical_data[MSP])
        else:
            raise Exception("Something went wrong!!")
        input_data[MSP] = list(map(lambda x: "Lower" if (clinical_data.loc[step00.get_patient(x), MSP] < threshold) else "Higher", list(input_data.index)))

    output_data = pandas.DataFrame(index=input_data.index)
    output_data["Stage"] = input_data["Stage"]
    for MSP in tqdm.tqdm(step00.sharing_columns):
        output_data[MSP] = input_data[MSP]
    with multiprocessing.Pool(args.cpus) as pool:
        output_data["L1"] = pool.map(L1, list(output_data.index))
        output_data["L2"] = pool.map(L2, list(output_data.index))
    print(output_data)

    figures = list()
    for distance, MSP in tqdm.tqdm(list(itertools.product(["L1", "L2"], step00.sharing_columns))):
        stage_list = list(filter(lambda x: (x in set(output_data.loc[(output_data[MSP] == "Lower"), "Stage"])) and (x in set(output_data.loc[(output_data[MSP] == "Higher"), "Stage"])), step00.long_sample_type_list[:-1]))

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        if len(stage_list) > 1:
            seaborn.violinplot(data=output_data, x="Stage", y=distance, order=stage_list, hue=MSP, hue_order=["Lower", "Higher"], cut=1, linewidth=5, ax=ax)
            statannotations.Annotator.Annotator(ax, list(map(lambda x: ((x, "Lower"), (x, "Higher")), stage_list)), data=output_data, x="Stage", y=distance, order=stage_list, hue=MSP, hue_order=["Lower", "Higher"]).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()
        else:
            seaborn.violinplot(data=output_data.loc[(output_data["Stage"] == stage_list[0])], x=MSP, y=distance, order=["Lower", "Higher"], palette={"Lower": "tab:blue", "Higher": "tab:red"}, cut=1, linewidth=5, ax=ax)
            statannotations.Annotator.Annotator(ax, [("Lower", "Higher")], data=output_data.loc[(output_data["Stage"] == stage_list[0])], x=MSP, y=distance, order=["Lower", "Higher"]).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

        matplotlib.pyplot.ylabel(f"{distance} distance")
        matplotlib.pyplot.tight_layout()

        figures.append(f"{distance}_{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
