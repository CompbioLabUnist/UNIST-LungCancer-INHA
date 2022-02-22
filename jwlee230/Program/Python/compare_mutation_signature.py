"""
compare_mutation_signature.py: compare mutation intersection proportion vs. mutation signature
"""
import argparse
import itertools
import tarfile
import multiprocessing
import matplotlib
import matplotlib.pyplot
import scipy.stats
import seaborn
import pandas
import tqdm
import step00


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

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(sorted(patients))

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))
    sample_list = list(map(step00.get_id, args.input))
    precancer_sample_list = list(filter(lambda x: step00.get_long_sample_type(x) != "Primary", sample_list))
    print(sample_list)

    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
        mutect_data = mutect_data.loc[(mutect_data["Variant_Classification"].isin(step00.nonsynonymous_mutations))]
        mutect_data["Tumor_Sample_Barcode"] = list(map(lambda x: x.split(".")[0], mutect_data["Tumor_Sample_Barcode"]))
        mutect_data["Stage"] = pool.map(step00.get_long_sample_type, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Patient"] = pool.map(step00.get_patient, mutect_data["Tumor_Sample_Barcode"])
    print(mutect_data)

    signature_data = pandas.read_csv(args.signature, sep="\t", index_col="Samples")
    for index in tqdm.tqdm(list(signature_data.index)):
        signature_data.loc[index, :] = signature_data.loc[index, :] / sum(signature_data.loc[index, :])
    print(signature_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    output_data = pandas.DataFrame(index=precancer_sample_list)
    output_data["Stage"] = list(map(step00.get_long_sample_type, list(output_data.index)))

    output_data["Proportion"] = 0
    for sample in tqdm.tqdm(precancer_sample_list):
        precancer_set = set(mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == sample), ["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]].itertuples(index=False, name=None))
        primary_set = set(mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == step00.get_paired_primary(sample)), ["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]].itertuples(index=False, name=None))

        output_data.loc[sample, "Proportion"] = len(precancer_set & primary_set) / len(precancer_set | primary_set)
    print(output_data)

    for signature in tqdm.tqdm(list(signature_data.columns)):
        output_data[signature + "_Precancer"] = signature_data.loc[precancer_sample_list, signature]
        output_data[signature + "_Primary"] = list(map(lambda x: signature_data.loc[step00.get_paired_primary(x), signature], precancer_sample_list))
    print(output_data)

    signature_list = list(output_data.columns)[2:]

    figures = list()

    for stage, signature in tqdm.tqdm(list(itertools.product(step00.long_sample_type_list, signature_list))):
        tmp_data = output_data.loc[(output_data["Stage"] == stage), ["Proportion", signature]]
        if tmp_data.shape[0] < 2:
            continue

        r, p = scipy.stats.pearsonr(tmp_data[signature], tmp_data["Proportion"])

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        g = seaborn.jointplot(data=tmp_data, x=signature, y="Proportion", kind="reg", height=24, ratio=6, xlim=(-0.1, 1.1), ylim=(-0.1, 1.1))
        g.fig.text(0.5, 0.75, "r={0:.3f}, p={1:.3f}".format(r, p), color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")
        g.set_axis_labels("{0} in {1}".format(signature.split("_")[0], signature.split("_")[1]), "Shared Proportion of non-SYN mutations")

        figures.append("{1}-{0}.pdf".format(signature, stage))
        g.savefig(figures[-1])
        matplotlib.pyplot.close()

    for signature in tqdm.tqdm(signature_list):
        tmp_data = output_data.loc[:, ["Proportion", signature]]

        if tmp_data.shape[0] < 2:
            continue

        r, p = scipy.stats.pearsonr(tmp_data[signature], tmp_data["Proportion"])

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        g = seaborn.jointplot(data=tmp_data, x=signature, y="Proportion", kind="reg", height=24, ratio=6, xlim=(-0.1, 1.1), ylim=(-0.1, 1.1))
        g.fig.text(0.5, 0.75, "r={0:.3f}, p={1:.3f}".format(r, p), color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")
        g.set_axis_labels("{0} in {1}".format(signature.split("_")[0], signature.split("_")[1]), "Shared Proportion of non-SYN mutations")

        figures.append("Precancer-{0}.pdf".format(signature))
        g.savefig(figures[-1])
        matplotlib.pyplot.close()

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
