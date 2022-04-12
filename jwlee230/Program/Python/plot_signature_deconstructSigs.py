"""
plot_signature_deconstructSigs.py: violin plot cancer signature from deconstructSigs by stage
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


def draw_violin(signature: str) -> pandas.DataFrame:
    fig, ax = matplotlib.pyplot.subplots(figsize=(7 * len(order), 24))

    seaborn.violinplot(data=input_data, x="Stage", y=signature, order=order, inner="box", ax=ax)
    statannotations.Annotator.Annotator(ax, list(itertools.combinations(order, 2)), data=input_data, x="Stage", y=signature, order=order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.ylabel("Count")
    matplotlib.pyplot.title(signature)

    fig_name = "{0}.pdf".format(signature)
    fig.savefig(os.path.join(step00.tmpfs, fig_name))
    matplotlib.pyplot.close(fig)

    outputs = [signature]
    for control, case in list(itertools.combinations(order, 2)):
        outputs.append(scipy.stats.mannwhitneyu(list(input_data.loc[(input_data["Stage"] == control), signature]), list(input_data.loc[(input_data["Stage"] == case), signature]))[1])

    return pandas.DataFrame(data=outputs, index=["Signature"] + list(map(lambda x: x[0] + "-" + x[1], list(itertools.combinations(order, 2))))).T


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Signature TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical data must end with .csv!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data.columns = list(map(lambda x: x.replace("Signature.", "SBS"), list(input_data.columns)))
    signatures = list(filter(lambda x: len(set(input_data[x])) > 1, list(input_data.columns)))
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

    input_data = input_data.loc[sorted(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index)), key=step00.sorting_by_type), :]
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    for stage in tqdm.tqdm(set(input_data["Stage"])):
        if len(input_data.loc[(input_data["Stage"] == stage)]) < 3:
            input_data = input_data.loc[~(input_data["Stage"] == stage)]
    order = list(filter(lambda x: x in set(input_data["Stage"]), step00.long_sample_type_list))
    print(input_data)
    print(order)

    with multiprocessing.Pool(args.cpus) as pool:
        output_data = pandas.concat(objs=pool.map(draw_violin, signatures), join="outer", ignore_index=True, axis="index").set_index(keys="Signature", verify_integrity=True)
    output_data.to_csv(os.path.join(step00.tmpfs, "output.tsv"), sep="\t", float_format="{:.2e}".format)
    print(output_data)

    with tarfile.open(name=args.output, mode="w") as tar:
        tar.add(os.path.join(step00.tmpfs, "output.tsv"), arcname="output.tsv")
        for file in tqdm.tqdm(list(output_data.index)):
            tar.add(os.path.join(step00.tmpfs, file + ".pdf"), arcname=file + ".pdf")
