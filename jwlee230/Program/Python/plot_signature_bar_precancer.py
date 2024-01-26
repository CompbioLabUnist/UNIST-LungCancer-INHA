"""
plot_signature_bar.py: Plot signature in bar graph comparing precancer vs. primary
"""
import argparse
import itertools
import tarfile
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import statannotations.Annotator
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Signature TSV file (not necessarily TSV)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--percentage", help="Percentage of patients to include", type=float, default=0.1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_type = parser.add_mutually_exclusive_group(required=True)
    group_type.add_argument("--absolute", help="Absolute bar graph", action="store_true", default=False)
    group_type.add_argument("--relative", help="Relative bar graph", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical data must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif not (0.0 < args.percentage < 0.5):
        raise ValueError("Percentage must be (0.0, 0.5)!!")

    input_data = pandas.DataFrame()
    for input_file in tqdm.tqdm(args.input):
        d = pandas.read_csv(input_file, sep="\t", index_col=0)
        if args.relative:
            for index in list(d.index):
                d.loc[index, :] = d.loc[index, :] / sum(d.loc[index, :])
        input_data = input_data.join(d, how="outer")
    signatures = list(input_data.columns)
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
    print(sorted(patients))

    input_data = input_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index))), :]
    print(input_data)

    samples = list(input_data.index)
    precancer_samples = list(filter(lambda x: (step00.get_paired_primary(x) in samples) and (step00.get_long_sample_type(x) != "Primary"), samples))
    print(precancer_samples)

    output_data = pandas.DataFrame([(signature, sample, input_data.loc[sample, signature]) for signature, sample in tqdm.tqdm(list(itertools.product(signatures, samples)))], columns=["Signature", "Sample", "Value"])
    output_data["PRE/PRI"] = list(map(lambda x: "Primary" if step00.get_long_sample_type(x) == "Primary" else "Precancer", output_data["Sample"]))
    print(output_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    for MSP, signature in tqdm.tqdm(list(itertools.product(step00.sharing_columns, signatures))):
        lower_bound, higher_bound = numpy.quantile(clinical_data[MSP], args.percentage), numpy.quantile(clinical_data[MSP], 1.0 - args.percentage)

        lower_precancer_list = list(clinical_data.loc[(clinical_data[MSP] <= lower_bound), f"{MSP}-sample"])
        higher_precancer_list = list(clinical_data.loc[(clinical_data[MSP] >= higher_bound), f"{MSP}-sample"])

        lower_primary_list = list(map(step00.get_paired_primary, lower_precancer_list))
        higher_primary_list = list(map(step00.get_paired_primary, higher_precancer_list))

        fig, axs = matplotlib.pyplot.subplots(figsize=(18, 36), nrows=2)

        seaborn.violinplot(data=output_data.loc[(output_data["Sample"].isin(lower_precancer_list + lower_primary_list))], x="Signature", y="Value", hue="PRE/PRI", order=[signature], hue_order=["Precancer", "Primary"], palette={"Precancer": "tab:pink", "Primary": "gray"}, inner="box", linewidth=5, cut=1, ax=axs[0])
        statannotations.Annotator.Annotator(axs[0], [((signature, "Precancer"), (signature, "Primary"))], data=output_data.loc[(output_data["Sample"].isin(lower_precancer_list + lower_primary_list))], x="Signature", y="Value", hue="PRE/PRI", order=[signature], hue_order=["Precancer", "Primary"]).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()

        seaborn.violinplot(data=output_data.loc[(output_data["Sample"].isin(higher_precancer_list + higher_primary_list))], x="Signature", y="Value", hue="PRE/PRI", order=[signature], hue_order=["Precancer", "Primary"], palette={"Precancer": "tab:pink", "Primary": "gray"}, inner="box", linewidth=5, cut=1, ax=axs[1])
        statannotations.Annotator.Annotator(axs[1], [((signature, "Precancer"), (signature, "Primary"))], data=output_data.loc[(output_data["Sample"].isin(higher_precancer_list + higher_primary_list))], x="Signature", y="Value", hue="PRE/PRI", order=[signature], hue_order=["Precancer", "Primary"]).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()

        axs[0].set_xlabel("")
        axs[0].set_ylabel(f"{signature} in MSP-lower")
        axs[1].set_xlabel("")
        axs[1].set_ylabel(f"{signature} in MSP-higher")

        matplotlib.pyplot.tight_layout()

        figures.append(f"{MSP}_{signature}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure)
