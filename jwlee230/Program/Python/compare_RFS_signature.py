"""
compare_RFS_signature.py: compare RFS (Recurrence Free Survival) vs. mutational signature
"""
import argparse
import itertools
import tarfile
import matplotlib
import matplotlib.pyplot
import scipy.stats
import seaborn
import pandas
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("signature", help="Mutation signature TSV file (not necessarily TSV)", type=str)
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(sorted(patients))

    signature_data = pandas.read_csv(args.signature, sep="\t", index_col="Samples")
    for index in tqdm.tqdm(list(signature_data.index)):
        signature_data.loc[index, :] = signature_data.loc[index, :] / sum(signature_data.loc[index, :])
    print(signature_data)

    sample_list = list(filter(lambda x: step00.get_patient(x) in patients, list(signature_data.index)))
    print(sample_list)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    output_data = pandas.DataFrame(index=sample_list)
    output_data["Stage"] = list(map(step00.get_long_sample_type, list(output_data.index)))
    output_data["Recurrence-Free Survival"] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), "Recurrence-Free Survial"], list(output_data.index)))
    print(output_data)

    signature_list = list(signature_data.columns)

    for signature in tqdm.tqdm(signature_list):
        output_data[signature] = signature_data.loc[sample_list, signature]
    print(output_data)

    figures = list()

    for stage, signature in tqdm.tqdm(list(itertools.product(step00.long_sample_type_list, signature_list))):
        tmp_data = output_data.loc[(output_data["Stage"] == stage), ["Recurrence-Free Survival", signature]]
        if tmp_data.shape[0] < 3:
            continue

        r, p = scipy.stats.pearsonr(tmp_data[signature], tmp_data["Recurrence-Free Survival"])

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        g = seaborn.jointplot(data=tmp_data, x=signature, y="Recurrence-Free Survival", kind="reg", height=24, ratio=6, xlim=(-0.1, 1.1))
        g.fig.text(0.5, 0.75, "r={0:.3f}, p={1:.3f}".format(r, p), color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")
        g.plot_marginals(seaborn.histplot, kde=True, stat="probability", multiple="stack")
        g.set_axis_labels("{0} proportion in {1}".format(signature, stage), "Recurrence-Free Survival (days)")

        figures.append("{1}-{0}.pdf".format(signature, stage))
        g.savefig(figures[-1])
        matplotlib.pyplot.close()

    for signature in tqdm.tqdm(signature_list):
        tmp_data = output_data.loc[:, ["Recurrence-Free Survival", signature]]

        if tmp_data.shape[0] < 3:
            continue

        r, p = scipy.stats.pearsonr(tmp_data[signature], tmp_data["Recurrence-Free Survival"])

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        g = seaborn.jointplot(data=tmp_data, x=signature, y="Recurrence-Free Survival", kind="reg", height=24, ratio=6, xlim=(-0.1, 1.1))
        g.fig.text(0.5, 0.75, "r={0:.3f}, p={1:.3f}".format(r, p), color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")
        g.plot_marginals(seaborn.histplot, kde=True, stat="probability", multiple="stack")
        g.set_axis_labels("{0} proportion".format(signature), "Recurrence-Free Survival (days)")

        figures.append("All-{0}.pdf".format(signature))
        g.savefig(figures[-1])
        matplotlib.pyplot.close()

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
