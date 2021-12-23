"""
draw_violin_plots_clinical.py: draw violin plots upon DEG with clinical data
"""
import argparse
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import statannotations.Annotator
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TPM TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data data CSV file", type=str)
    parser.add_argument("cgc", help="CGC CSV files", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--compare", help="Comparison grouping (type, control, case)", type=str, nargs=3, default=["Recurrence", "NO", "YES"])
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .tsv!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.cgc.endswith(".csv"):
        raise ValueError("CGC must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        control_patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC") & (clinical_data[args.compare[0]] == args.compare[1])].index)
        case_patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC") & (clinical_data[args.compare[0]] == args.compare[2])].index)
    elif args.ADC:
        control_patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC") & (clinical_data[args.compare[0]] == args.compare[1])].index)
        case_patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC") & (clinical_data[args.compare[0]] == args.compare[2])].index)
    else:
        raise Exception("Something went wrong!!")
    patients = control_patients | case_patients
    print(sorted(control_patients))
    print(sorted(case_patients))

    cgc_data = pandas.read_csv(args.cgc, index_col="Gene Symbol")
    gene_set = set(cgc_data.index)
    print(cgc_data)

    ylabel = args.input.split("/")[-1].split(".")[1]

    input_data = pandas.read_csv(args.input, sep="\t", index_col="gene_name").T

    samples = list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index)))
    gene_set &= set(input_data.columns)

    input_data = input_data.loc[samples, sorted(gene_set)]
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    input_data[args.compare[0]] = list(map(lambda x: args.compare[1] if step00.get_patient(x) in control_patients else args.compare[2], list(input_data.index)))
    print(input_data)

    stage_order = list(filter(lambda x: (x in set(map(step00.get_long_sample_type, input_data.loc[(input_data[args.compare[0]] == args.compare[1])].index))) and (x in set(map(step00.get_long_sample_type, input_data.loc[(input_data[args.compare[0]] == args.compare[2])].index))), step00.long_sample_type_list))
    print(stage_order)

    figures = list()

    for gene in tqdm.tqdm(gene_set):
        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        seaborn.violinplot(data=input_data, x="Stage", y=gene, order=stage_order, ax=ax, hue=args.compare[0], hue_order=args.compare[1:])
        statannotations.Annotator.Annotator(ax, [((stage, args.compare[1]), (stage, args.compare[2])) for stage in stage_order], data=input_data, x="Stage", y=gene, order=stage_order, hue=args.compare[0], hue_order=args.compare[1:]).configure(test="Mann-Whitney", text_format="star", loc="inside", verbose=0).apply_and_annotate()

        matplotlib.pyplot.ylabel(ylabel)
        matplotlib.pyplot.title(gene)
        matplotlib.pyplot.tight_layout()

        figures.append(gene + ".pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(figures):
            tar.add(f, arcname=f)
