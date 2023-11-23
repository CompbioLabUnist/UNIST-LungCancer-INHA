"""
draw_bar_plot_deconvolution_precancer.py: draw bar plot from deconvolution results comparing precancer vs. primary
"""
import argparse
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Deconvolution result TSV file (not necessarily TSV)", type=str)
    parser.add_argument("clinical", help="Clinical data w/ Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    else:
        raise Exception("Something went wrong!!")
    print(clinical_data)

    patients = set(clinical_data.index)
    print(patients)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data = input_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index))), :]
    print(input_data)

    samples = list(input_data.index)
    precancer_samples = list(filter(lambda x: (step00.get_long_sample_type(x) not in {"Normal", "Primary"}) and (step00.get_paired_primary(x) in samples), samples))
    print(len(precancer_samples))

    cells = sorted(input_data.columns, key=lambda x: sum(input_data[x]), reverse=True)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    figures = list()
    for precancer_sample in tqdm.tqdm(precancer_samples):
        fig, ax = matplotlib.pyplot.subplots(figsize=(64, 18))

        matplotlib.pyplot.bar(range(len(cells)), input_data.loc[precancer_sample, cells], width=-0.4, align="edge", color=step00.stage_color_code[step00.get_long_sample_type(precancer_sample)], label=step00.get_long_sample_type(precancer_sample))
        matplotlib.pyplot.bar(range(len(cells)), input_data.loc[step00.get_paired_primary(precancer_sample), cells], width=0.4, align="edge", color=step00.stage_color_code["Primary"], label="Primary")

        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.legend(loc="upper right")
        matplotlib.pyplot.xticks(range(len(cells)), cells, rotation="vertical")
        matplotlib.pyplot.ylabel("Cell proportion")
        matplotlib.pyplot.title(f"{step00.get_patient(precancer_sample)}")
        matplotlib.pyplot.tight_layout()

        figures.append(f"{precancer_sample}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure)
