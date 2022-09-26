"""
draw_violin_beta_PathSeq.py: draw violin plot with Beta-diversity of PathSeq
"""
import argparse
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import statannotations.Annotator
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PathSeq beta-diversity results TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs", type=int, default=1)

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
        raise ValueError("CPUS must be positive!!")

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
    samples = list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index)))
    input_data = input_data.loc[samples, samples]
    print(input_data)

    input_data["Patient"] = list(map(step00.get_patient, list(input_data.index)))
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    print(input_data)

    output_data = pandas.DataFrame(index=list(filter(lambda x: step00.get_long_sample_type(x) not in {"Primary"}, samples)))
    output_data["Beta-diversity"] = list(map(lambda x: input_data.loc[x, list(input_data.loc[(input_data["Patient"] == step00.get_patient(x)) & (input_data["Stage"] == "Primary")].index)[0]], list(output_data.index)))
    output_data["Patient"] = list(map(step00.get_patient, list(output_data.index)))
    output_data["Stage"] = list(map(step00.get_long_sample_type, list(output_data.index)))
    print(output_data)

    for MSP in tqdm.tqdm(step00.sharing_columns):
        if args.median:
            threshold = numpy.median(clinical_data[MSP])
        elif args.mean:
            threshold = numpy.mean(clinical_data[MSP])
        else:
            raise ValueError("Something went wrong!!")

        output_data[MSP] = list(map(lambda x: "Lower" if (clinical_data.loc[step00.get_patient(x), MSP] < threshold) else "Higher", list(output_data.index)))
    print(output_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()

    for MSP in tqdm.tqdm(step00.sharing_columns):
        stage_list = list(filter(lambda x: (x in list(output_data.loc[(output_data[MSP] == "Lower"), "Stage"])) and (x in list(output_data.loc[(output_data[MSP] == "Higher"), "Stage"])), step00.long_sample_type_list))

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.violinplot(data=output_data, x="Stage", order=stage_list, y="Beta-diversity", hue=MSP, hue_order=["Lower", "Higher"], palette={"Lower": "tab:cyan", "Higher": "tab:red"}, inner="box", cut=1, ax=ax)
        statannotations.Annotator.Annotator(ax, [((stage, "Lower"), (stage, "Higher")) for stage in stage_list], data=output_data, x="Stage", order=stage_list, y="Beta-diversity", hue=MSP, hue_order=["Lower", "Higher"]).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

        matplotlib.pyplot.tight_layout()

        figures.append(f"{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
