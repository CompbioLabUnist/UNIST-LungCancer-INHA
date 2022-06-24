"""
plot_signature_2_relative.py: violin plot cancer signature with relative by stage with clinical data
"""
import argparse
import itertools
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
compare_order: typing.List[typing.Tuple[typing.Tuple[str, str], typing.Tuple[str, str]]] = list()
hue_order: typing.List[str] = list()


def draw_violin(signature: str, clinical: str) -> pandas.DataFrame:
    try:
        stat, p = scipy.stats.kruskal(*[input_data.loc[(input_data["Stage"] == stage) & (input_data[clinical] == clinical_value), signature] for stage, clinical_value in itertools.product(order, hue_order)])
    except ValueError:
        _, p = 0.0, 1.0

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.violinplot(data=input_data, x="Stage", y=signature, order=order, hue=clinical, hue_order=hue_order, inner="box", cut=1, ax=ax)
    statannotations.Annotator.Annotator(ax, compare_order, data=input_data, x="Stage", y=signature, order=order, hue=clinical, hue_order=hue_order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.ylabel("Proportion")
    matplotlib.pyplot.title(f"{signature}: Kruskal-Wallis p={p:.3f}")
    matplotlib.pyplot.tight_layout()

    fig_name = "{0}.pdf".format(signature)
    fig.savefig(os.path.join(step00.tmpfs, fig_name))
    matplotlib.pyplot.close(fig)

    outputs = [signature]
    for (s1, c1), (s2, c2) in compare_order:
        outputs.append(scipy.stats.mannwhitneyu(list(input_data.loc[(input_data["Stage"] == s1) & (input_data[clinical] == c1), signature]), list(input_data.loc[(input_data["Stage"] == s2) & (input_data[clinical] == c2), signature]))[1])

    return pandas.DataFrame(data=outputs, index=["Signature"] + ["{0}: {1}-{2}".format(s1, c1, c2) for (s1, c1), (s2, c2) in compare_order]).T


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Signature TSV file (not necessarily TSV)", type=str)
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--compare", help="Comparison grouping (type, control, case)", type=str, nargs="+", default=["Recurrence", "NO", "YES"])
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".csv"):
        raise ValueError("Clinical data must end with .csv!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = pandas.read_csv(args.input, sep="\t", index_col="Samples")
    signatures = list(input_data.columns)
    print(input_data)

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC") & (clinical_data[args.compare[0]].isin(args.compare[1:]))].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC") & (clinical_data[args.compare[0]].isin(args.compare[1:]))].index)
    else:
        raise Exception("Something went wrong!!")

    input_data = input_data.loc[sorted(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index)), key=step00.sorting_by_type), :]
    input_data["Total"] = input_data.sum(axis="columns")
    for index in list(input_data.index):
        input_data.loc[index, :] = input_data.loc[index, :] / input_data.loc[index, "Total"]
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    input_data[args.compare[0]] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), args.compare[0]], list(input_data.index)))
    print(input_data)

    for stage in tqdm.tqdm(set(input_data["Stage"])):
        flag = False
        for clinical in args.compare[1:]:
            if input_data.loc[(input_data["Stage"] == stage) & (input_data[args.compare[0]] == clinical)].empty:
                flag = True
                break

        if flag:
            input_data = input_data.loc[~(input_data["Stage"] == stage)]
    hue_order = args.compare[1:]
    order = list(filter(lambda x: x in list(input_data["Stage"]), step00.long_sample_type_list))
    compare_order = list(filter(lambda x: not input_data.loc[(input_data["Stage"] == x[0][0]) & (input_data[args.compare[0]] == x[0][1])].empty and not input_data.loc[(input_data["Stage"] == x[1][0]) & (input_data[args.compare[0]] == x[1][1])].empty, [((s, c1), (s, c2)) for c1, c2 in itertools.combinations(hue_order, 2) for s in order]))
    print(order)
    print(hue_order)
    print(compare_order)

    with multiprocessing.Pool(args.cpus) as pool:
        output_data = pandas.concat(objs=pool.starmap(draw_violin, [(signature, args.compare[0]) for signature in signatures]), join="outer", ignore_index=True, axis="index").set_index(keys="Signature", verify_integrity=True)
    output_data.to_csv(os.path.join(step00.tmpfs, "output.tsv"), sep="\t", float_format="{:.2e}".format)
    print(output_data)

    with tarfile.open(name=args.output, mode="w") as tar:
        tar.add(os.path.join(step00.tmpfs, "output.tsv"), arcname="output.tsv")
        for file in tqdm.tqdm(list(output_data.index)):
            tar.add(os.path.join(step00.tmpfs, file + ".pdf"), arcname=file + ".pdf")
