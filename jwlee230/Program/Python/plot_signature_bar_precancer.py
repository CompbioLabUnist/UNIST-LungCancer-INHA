"""
plot_signature_bar.py: Plot signature in bar graph comparing precancer vs. primary
"""
import argparse
import tarfile
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import pandas
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Signature TSV file (not necessarily TSV)", type=str)
    parser.add_argument("clinical", help="Clinical data data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_type = parser.add_mutually_exclusive_group(required=True)
    group_type.add_argument("--absolute", help="Absolute bar graph", action="store_true", default=False)
    group_type.add_argument("--relative", help="Relative bar graph", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".csv"):
        raise ValueError("Clinical data must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col="Samples")
    signatures = list(input_data.columns)
    print(input_data)

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(sorted(patients))

    input_data = input_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index))), :]
    signatures = list(input_data.columns)
    input_data["Total"] = input_data.sum(axis="columns")
    print(input_data)
    print(signatures)

    if args.absolute:
        input_data.sort_values(by=sorted(signatures, key=lambda x: sum(input_data.loc[:, x]), reverse=True), ascending=False, inplace=True)
        input_data.sort_values(by="Total", kind="stable", ascending=False, inplace=True)
        input_data = input_data.loc[:, signatures]
    elif args.relative:
        for index in tqdm.tqdm(list(input_data.index)):
            input_data.loc[index, :] = input_data.loc[index, :] / input_data.loc[index, "Total"]
        input_data.sort_values(by=sorted(signatures, key=lambda x: sum(input_data.loc[:, x]), reverse=True), ascending=False, inplace=True)
        input_data = input_data.loc[:, signatures]
    else:
        raise Exception("Something went wrong!!")
    print(input_data)

    samples = list(input_data.index)
    precancer_samples = list(filter(lambda x: (step00.get_paired_primary(x) in samples) and (step00.get_long_sample_type(x) != "Primary"), samples))
    print(precancer_samples)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    figures = list()
    for precancer_sample in tqdm.tqdm(precancer_samples):
        fig, ax = matplotlib.pyplot.subplots(figsize=(64, 18))

        matplotlib.pyplot.bar(range(len(signatures)), list(input_data.loc[precancer_sample, signatures]), align="edge", width=-0.4, color=step00.get_color_by_type(precancer_sample), label=step00.get_long_sample_type(precancer_sample))
        matplotlib.pyplot.bar(range(len(signatures)), list(input_data.loc[step00.get_paired_primary(precancer_sample), signatures]), align="edge", width=0.4, color=step00.stage_color_code["Primary"], label="Primary")

        if args.absolute:
            matplotlib.pyplot.ylabel("Counts")
        elif args.relative:
            matplotlib.pyplot.ylabel("Proportion")
        else:
            raise Exception("Something went wrong!!")

        matplotlib.pyplot.xlabel(f"Signature in {step00.get_patient(precancer_sample)}")
        matplotlib.pyplot.xticks(range(len(signatures)), signatures, rotation="vertical")
        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.legend()
        matplotlib.pyplot.tight_layout()

        figures.append(f"{precancer_sample}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure)
