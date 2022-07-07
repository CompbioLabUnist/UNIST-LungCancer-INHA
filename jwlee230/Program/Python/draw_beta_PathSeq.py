"""
draw_beta_PathSeq.py: draw beta-diversity by stage
"""
import argparse
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import skbio.stats.distance
import sklearn.manifold
import sklearn.preprocessing
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PathSeq beta-diversity results TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUS must be positive!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    samples = list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index)))
    input_data = input_data.loc[samples, samples]
    print(input_data)

    p = skbio.stats.distance.permanova(skbio.stats.distance.DistanceMatrix(numpy.ascontiguousarray(input_data.to_numpy()), ids=samples), grouping=list(map(step00.get_long_sample_type, samples)))["p-value"]
    print(p)

    output_data = pandas.DataFrame(sklearn.manifold.TSNE(init="pca", verbose=1, random_state=42, method="exact", n_jobs=args.cpus).fit_transform(input_data), index=samples, columns=["tSNE1", "tSNE2"])
    for column in tqdm.tqdm(list(output_data.columns)):
        output_data[column] = sklearn.preprocessing.scale(output_data[column])
    output_data["Subtype"] = list(map(step00.get_long_sample_type, samples))
    print(output_data)

    order = list(filter(lambda x: x in set(output_data["Subtype"]), step00.long_sample_type_list))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.scatterplot(data=output_data, x="tSNE1", y="tSNE2", hue="Subtype", style="Subtype", hue_order=order, palette=step00.stage_color_code, s=1000)

    matplotlib.pyplot.title(f"PERMANOVA p={p:.3f}")
    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
