"""
plot_expression_clinical_2.py: Plot gene expression TPM with cancer subtype with clinical data
"""
import argparse
import multiprocessing
import os.path
import tarfile
import typing
import matplotlib
import matplotlib.pyplot
import pandas
import scipy.stats
import seaborn
import statannotations.Annotator
import tqdm
import step00

input_data = pandas.DataFrame()
order: typing.List[str] = list()
compare_order: typing.List[str] = list()
hue_order: typing.List[str] = list()


def draw_violin(gene: str, clinical: str) -> pandas.DataFrame:
    fig, ax = matplotlib.pyplot.subplots(figsize=(7 * len(order), 24))

    seaborn.violinplot(data=input_data, x="Subtype", y=gene, order=order, hue=clinical, hue_order=hue_order, inner="box", ax=ax)
    statannotations.Annotator.Annotator(ax, list(map(lambda x: [(x, hue_order[0]), (x, hue_order[1])], compare_order)), data=input_data, x="Subtype", y=gene, order=order, hue=clinical, hue_order=hue_order).configure(test="Mann-Whitney", text_format="star", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.ylabel("TPM")
    matplotlib.pyplot.title(gene)

    fig_name = "{0}.pdf".format(gene)
    fig.savefig(os.path.join(step00.tmpfs, fig_name))
    matplotlib.pyplot.close(fig)

    outputs = [gene]
    for subtype in compare_order:
        outputs.append(scipy.stats.mannwhitneyu(list(input_data.loc[(input_data["Subtype"] == subtype) & (input_data[clinical] == hue_order[0]), gene]), list(input_data.loc[(input_data["Subtype"] == subtype) & (input_data[clinical] == hue_order[1]), gene]))[1])

    return pandas.DataFrame(data=outputs, index=["Gene"] + compare_order).T


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="RSEM TPM files", type=str)
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--compare", help="Comparison grouping (type, control, case)", type=str, nargs=3, default=["Recurrence", "NO", "YES"])
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .tsv!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical data must end with .csv!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = pandas.read_csv(args.input, sep="\t", index_col="gene_name")
    genes = list(input_data.index)
    print(input_data)

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
    patients = set(list(control_patients) + list(case_patients))
    print(sorted(control_patients))
    print(sorted(case_patients))

    input_data = input_data.loc[:, sorted(filter(lambda x: step00.get_patient(x) in patients, list(input_data.columns)), key=step00.sorting_by_type)].T
    input_data["Subtype"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    input_data[args.compare[0]] = list(map(lambda x: args.compare[1] if (step00.get_patient(x) in control_patients) else args.compare[2], list(input_data.index)))
    print(input_data)

    order = list(filter(lambda x: x in set(input_data["Subtype"]), step00.long_sample_type_list))
    compare_order = list(filter(lambda x: not input_data.loc[(input_data["Subtype"] == x) & (input_data[args.compare[0]] == args.compare[1])].empty and not input_data.loc[(input_data["Subtype"] == x) & (input_data[args.compare[0]] == args.compare[2])].empty, order))
    hue_order = args.compare[1:]
    print(order)
    print(compare_order)
    print(hue_order)

    with multiprocessing.Pool(args.cpus) as pool:
        output_data = pandas.concat(objs=pool.starmap(draw_violin, [(gene, args.compare[0]) for gene in genes]), join="outer", ignore_index=True, axis="index").set_index(keys="Gene", verify_integrity=True)
    output_data.to_csv(os.path.join(step00.tmpfs, "output.tsv"), sep="\t", float_format="{:.2e}".format)
    print(output_data)

    with tarfile.open(name=args.output, mode="w") as tar:
        tar.add(os.path.join(step00.tmpfs, "output.tsv"), arcname="output.tsv")
        for file in tqdm.tqdm(list(output_data.index)):
            tar.add(os.path.join(step00.tmpfs, file), arcname=file)
