"""
draw_tSNE_MutationSharedProportion_SigProfiler.py: draw t-SNE plot with Mutation Shared Proportion and Signature
"""
import argparse
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import sklearn.manifold
import sklearn.preprocessing
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("signature", help="Signature TXT (not necessarily TSV) file", type=str)
    parser.add_argument("clinical", help="Clinical data w/ Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_threshold = parser.add_mutually_exclusive_group(required=True)
    group_threshold.add_argument("--median", help="Use Median threshold", action="store_true", default=False)
    group_threshold.add_argument("--mean", help="Use Mean threshold", action="store_true", default=False)

    group_relative = parser.add_mutually_exclusive_group(required=True)
    group_relative.add_argument("--absolute", help="Use Absolute(Count) value", action="store_true", default=False)
    group_relative.add_argument("--relative", help="Use Relative(Proportion) value", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical data must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    print(clinical_data)

    input_data = pandas.read_csv(args.signature, sep="\t", index_col=0)
    if args.relative:
        for index in tqdm.tqdm(list(input_data.index)):
            input_data.loc[index, :] = input_data.loc[index, :] / sum(input_data.loc[index, :])
    print(input_data)

    signatures = list(input_data.columns)
    input_data["Patient"] = list(map(step00.get_patient, list(input_data.index)))
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    print(input_data)

    input_data = input_data.loc[(input_data["Patient"].isin(set(clinical_data.index)))]
    print(input_data)

    stage_list = list(filter(lambda x: x in set(input_data["Stage"]), step00.long_sample_type_list))
    print(stage_list)

    tsne_data = pandas.DataFrame(sklearn.manifold.TSNE(init="pca", verbose=1, random_state=0, method="exact", n_jobs=args.cpus).fit_transform(input_data.loc[:, signatures]), columns=["tSNE1", "tSNE2"], index=input_data.index)
    for column in tqdm.tqdm(list(tsne_data.columns)):
        tsne_data[column] = sklearn.preprocessing.scale(tsne_data[column])
    tsne_data["Patient"] = list(map(step00.get_patient, list(tsne_data.index)))
    tsne_data["Stage"] = list(map(step00.get_long_sample_type, list(tsne_data.index)))
    print(tsne_data)

    figures = list()
    for MSP in tqdm.tqdm(step00.sharing_columns):
        if args.median:
            threshold = numpy.median(clinical_data[MSP])
        elif args.mean:
            threshold = numpy.mean(clinical_data[MSP])
        else:
            raise Exception("Something went wrong!!")
        tsne_data[MSP] = list(map(lambda x: "Lower" if (clinical_data.loc[step00.get_patient(x), MSP] < threshold) else "Higher", list(tsne_data.index)))

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        for index in list(tsne_data.index):
            x1, y1 = tsne_data.loc[index, ["tSNE1", "tSNE2"]]
            x2, y2 = tsne_data.loc[step00.get_paired_primary(index), ["tSNE1", "tSNE2"]]

            matplotlib.pyplot.arrow(x1, y1, x2 - x1, y2 - y1, length_includes_head=True, color="tab:blue" if (tsne_data.loc[index, MSP] == "Lower") else "tab:red", alpha=0.3)

        seaborn.scatterplot(data=tsne_data, x="tSNE1", y="tSNE2", hue=MSP, style="Stage", palette={"Lower": "tab:blue", "Higher": "tab:red"}, style_order=stage_list, legend="full", s=1000, ax=ax)

        matplotlib.pyplot.title(f"Threshold: {threshold * 100:.2f}%")
        matplotlib.pyplot.tight_layout()

        figures.append(f"{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
