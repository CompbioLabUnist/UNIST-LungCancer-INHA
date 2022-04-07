"""
compare_SharedProportion_signature.py: compare mutation shared proportion vs. mutational signature
"""
import argparse
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import scipy.stats
import seaborn
import pandas
import tqdm
import step00

wanted_columns = ["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("signature", help="Mutation signature TSV file (not necessarily TSV)", type=str)
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    clinical_data: pandas.DataFrame = step00.get_clinical_data(args.clinical)
    clinical_data = clinical_data.loc[~(clinical_data["Volume_Doubling_Time"].isna())]
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(patients)

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))
    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
        mutect_data["Tumor_Sample_Barcode"] = pool.map(step00.get_id, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Patient"] = pool.map(step00.get_patient, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Stage"] = pool.map(step00.get_long_sample_type, mutect_data["Tumor_Sample_Barcode"])
    print(mutect_data)

    signature_data = pandas.read_csv(args.signature, sep="\t", index_col="Samples")
    signature_list = list(signature_data.columns)
    for index in tqdm.tqdm(list(signature_data.index)):
        signature_data.loc[index, :] = signature_data.loc[index, :] / sum(signature_data.loc[index, :])
    signature_data["Patient"] = list(map(step00.get_patient, list(signature_data.index)))
    signature_data["Stage"] = list(map(step00.get_long_sample_type, list(signature_data.index)))
    print(signature_data)

    patients &= set(mutect_data["Patient"])
    patients &= set(signature_data["Patient"])

    signature_data = signature_data.loc[signature_data["Patient"].isin(patients)]

    signature_data["Shared Proportion"] = 0.0
    for index in tqdm.tqdm(list(signature_data.index)):
        patient = step00.get_patient(index)
        stage = step00.get_long_sample_type(index)

        if stage == "Primary":
            signature_data.loc[index, "Shared Proportion"] = max(signature_data.loc[(signature_data["Patient"] == patient), "Shared Proportion"])
        else:
            patient_data = mutect_data.loc[mutect_data["Patient"] == patient]
            assert "Primary" in set(patient_data["Stage"])
            primary_set = set(patient_data.loc[patient_data["Stage"] == "Primary", wanted_columns].itertuples(index=False, name=None))

            precancer_set = set(patient_data.loc[patient_data["Stage"] == stage, wanted_columns].itertuples(index=False, name=None))
            proportion = len(primary_set & precancer_set) / len(primary_set)
            signature_data.loc[index, "Shared Proportion"] = proportion

    print(signature_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()

    for stage, signature in tqdm.tqdm(list(itertools.product(step00.long_sample_type_list, signature_list))):
        tmp_data = signature_data.loc[(signature_data["Stage"] == stage), ["Shared Proportion", signature]]
        if tmp_data.shape[0] < 3:
            continue

        r, p = scipy.stats.pearsonr(tmp_data[signature], tmp_data["Shared Proportion"])

        g = seaborn.jointplot(data=tmp_data, x=signature, y="Shared Proportion", kind="reg", height=24, ratio=6, xlim=(-0.1, 1.1), ylim=(-0.1, 1.1))
        g.fig.text(0.5, 0.75, "r={0:.3f}, p={1:.3f}".format(r, p), color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")
        g.plot_marginals(seaborn.histplot, kde=True, stat="probability", multiple="stack")
        g.set_axis_labels("{0} proportion in {1}".format(signature, stage), "Shared Proportion")

        figures.append("{1}-{0}.pdf".format(signature, stage))
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

    for signature in tqdm.tqdm(signature_list):
        tmp_data = signature_data.loc[:, ["Shared Proportion", signature]]

        if tmp_data.shape[0] < 3:
            continue

        r, p = scipy.stats.pearsonr(tmp_data[signature], tmp_data["Shared Proportion"])

        g = seaborn.jointplot(data=tmp_data, x=signature, y="Shared Proportion", kind="reg", height=24, ratio=6)
        g.fig.text(0.5, 0.75, "r={0:.3f}, p={1:.3f}".format(r, p), color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")
        g.plot_marginals(seaborn.histplot, kde=True, stat="probability", multiple="stack")
        g.set_axis_labels("{0} proportion".format(signature), "Shared Proportion")

        figures.append("All-{0}.pdf".format(signature))
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
