"""
draw_bar_plot_deconvolution_precancer.py: draw bar plot from deconvolution results comparing precancer vs. primary
"""
import argparse
import itertools
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import statannotations.Annotator
import tqdm
import step00

replaced = " ()/"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Deconvolution result TSV file (not necessarily TSV)", type=str)
    parser.add_argument("clinical", help="Clinical data w/ Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--percentage", help="Percentage of patients to include", type=float, default=0.25)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif not (0.0 < args.percentage < 0.5):
        raise ValueError("Percentage must be (0.0, 0.5)!!")

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
    cells = sorted(input_data.columns)
    print(input_data)

    samples = list(input_data.index)
    precancer_samples = list(filter(lambda x: (step00.get_long_sample_type(x) not in {"Normal", "Primary"}) and (step00.get_paired_primary(x) in samples), samples))
    print(len(precancer_samples))

    output_data = pandas.DataFrame([(cell, sample, input_data.loc[sample, cell]) for cell, sample in tqdm.tqdm(list(itertools.product(cells, samples)))], columns=["Cell", "Sample", "Proportion"])
    output_data["PRE/PRI"] = list(map(lambda x: "Primary" if step00.get_long_sample_type(x) == "Primary" else "Precancer", output_data["Sample"]))
    print(output_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    for MSP, cell in tqdm.tqdm(list(itertools.product(step00.sharing_columns, cells))):
        lower_bound, higher_bound = numpy.quantile(clinical_data[MSP], args.percentage), numpy.quantile(clinical_data[MSP], 1.0 - args.percentage)

        lower_precancer_list = list(clinical_data.loc[(clinical_data[MSP] < lower_bound), f"{MSP}-sample"])
        higher_precancer_list = list(clinical_data.loc[(clinical_data[MSP] >= higher_bound), f"{MSP}-sample"])

        lower_primary_list = list(map(step00.get_paired_primary, lower_precancer_list))
        higher_primary_list = list(map(step00.get_paired_primary, higher_precancer_list))

        lower_patients = list(map(step00.get_patient, lower_precancer_list))
        higher_patients = list(map(step00.get_patient, higher_precancer_list))

        drawing_data = output_data.loc[(output_data["Sample"].isin(lower_precancer_list + lower_primary_list + higher_precancer_list + higher_primary_list)) & (output_data["Cell"] == cell)].copy()
        drawing_data[MSP] = list(map(lambda x: "Lower" if (step00.get_patient(x) in lower_patients) else "Higher", drawing_data["Sample"]))

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.violinplot(data=drawing_data, x=MSP, y="Proportion", order=["Lower", "Higher"], hue="PRE/PRI", hue_order=["Precancer", "Primary"], palette={"Precancer": "tab:pink", "Primary": "gray"}, inner="box", linewidth=5, cut=1, ax=ax)
        statannotations.Annotator.Annotator(ax, [(("Lower", "Precancer"), ("Lower", "Primary")), (("Higher", "Precancer"), ("Higher", "Primary")), (("Lower", "Precancer"), ("Higher", "Precancer")), (("Lower", "Primary"), ("Higher", "Primary"))], data=drawing_data, x=MSP, y="Proportion", order=["Lower", "Higher"], hue="PRE/PRI", hue_order=["Precancer", "Primary"]).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()

        matplotlib.pyplot.ylabel(f"{cell} proportion")
        matplotlib.pyplot.tight_layout()

        cell_name = cell[:]
        for r in replaced:
            cell_name = cell_name.replace(r, "")

        figures.append(f"{MSP}-{cell_name}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure)
