"""
compare_mutation_shared_proportion.py: Compare mutation shared proportion
"""
import argparse
import collections
import multiprocessing
import matplotlib
import matplotlib.pyplot
import numpy
import seaborn
import pandas
import tqdm
import step00

wanted_columns = ["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]
watching = "Hugo_Symbol"


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--threshold", help="Threshold to use", type=int, default=20)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_strategy = parser.add_mutually_exclusive_group(required=True)
    group_strategy.add_argument("--median", help="Median division", action="store_true", default=False)
    group_strategy.add_argument("--mean", help="Mean division", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif args.threshold <= 5:
        raise ValueError("Threshold must be greater than 5!!")

    clinical_data: pandas.DataFrame = step00.get_clinical_data(args.clinical)
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

    patients &= set(mutect_data["Patient"])

    clinical_data["Shared Proportion"] = 0.0
    for patient in tqdm.tqdm(patients):
        patient_data = mutect_data.loc[mutect_data["Patient"] == patient]

        stage_set = set(patient_data["Stage"])
        assert "Primary" in stage_set
        primary_set = set(patient_data.loc[patient_data["Stage"] == "Primary", wanted_columns].itertuples(index=False, name=None))

        for stage in stage_set:
            if stage == "Primary":
                continue

            precancer_set = set(patient_data.loc[patient_data["Stage"] == stage, wanted_columns].itertuples(index=False, name=None))
            proportion = len(primary_set & precancer_set) / len(primary_set)

            clinical_data.loc[patient, "Shared Proportion"] = max(clinical_data.loc[patient, "Shared Proportion"], proportion)
    print(clinical_data)

    if args.median:
        cutting = numpy.median(clinical_data["Shared Proportion"])
    elif args.mean:
        cutting = numpy.mean(clinical_data["Shared Proportion"])
    else:
        raise Exception("Something went wrong!!")

    lower_data = clinical_data.loc[(clinical_data["Shared Proportion"] < cutting)]
    higher_data = clinical_data.loc[(clinical_data["Shared Proportion"] >= cutting)]

    lower_genes: collections.Counter = collections.Counter(mutect_data.loc[(mutect_data["Patient"].isin(set(lower_data.index))) & (mutect_data["Variant_Classification"].isin(step00.nonsynonymous_mutations))].drop_duplicates(subset=wanted_columns + ["Tumor_Sample_Barcode"]).loc[:, "Hugo_Symbol"])
    higher_genes: collections.Counter = collections.Counter(mutect_data.loc[(mutect_data["Patient"].isin(set(higher_data.index))) & (mutect_data["Variant_Classification"].isin(step00.nonsynonymous_mutations))].drop_duplicates(subset=wanted_columns + ["Tumor_Sample_Barcode"]).loc[:, "Hugo_Symbol"])
    total_genes = lower_genes + higher_genes

    lower_gene_names = set(lower_genes.keys())
    higher_gene_names = set(higher_genes.keys())

    output_data = pandas.DataFrame([(gene, count, "Lower") for gene, count in lower_genes.most_common()] + [(gene, count, "Higher") for gene, count in higher_genes.most_common()], columns=["Gene", "Count", "Lower/Higher"])
    print(output_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, axs = matplotlib.pyplot.subplots(figsize=(32, 48), nrows=3, sharey=True)

    genes = sorted(set(lower_gene_names) & set(higher_gene_names), key=lambda x: numpy.mean(output_data.loc[(output_data["Gene"] == x), "Count"]), reverse=True)[:args.threshold]
    tmp_data = output_data.loc[(output_data["Gene"].isin(genes))].sort_values(by="Count", ascending=False)
    seaborn.barplot(data=tmp_data, x="Gene", y="Count", order=genes, hue="Lower/Higher", ci=None, palette={"Lower": "tab:cyan", "Higher": "tab:red"}, ax=axs[0])
    axs[0].set_title(f"Both {len(clinical_data)} patients")
    axs[0].set_xticklabels(genes, rotation="vertical")
    axs[0].set_xlabel("")

    genes = sorted(set(lower_gene_names) - set(higher_gene_names))
    print(len(genes))
    tmp_data = output_data.loc[(output_data["Gene"].isin(genes)) & (output_data["Lower/Higher"] == "Lower")].sort_values(by="Count", ascending=False).iloc[:args.threshold, :]
    seaborn.barplot(data=tmp_data, x="Gene", y="Count", color="tab:cyan", ci=None, ax=axs[1])
    axs[1].set_title(f"Lower {len(lower_data)} patients")
    axs[1].set_xticklabels(tmp_data["Gene"], rotation="vertical")
    axs[1].set_xlabel("")

    genes = sorted(set(higher_gene_names) - set(lower_gene_names))
    print(len(genes))
    tmp_data = output_data.loc[(output_data["Gene"].isin(genes)) & (output_data["Lower/Higher"] == "Higher")].sort_values(by="Count", ascending=False).iloc[:args.threshold, :]
    seaborn.barplot(data=tmp_data, x="Gene", y="Count", color="tab:red", ci=None, ax=axs[2])
    axs[2].set_title(f"Higher {len(higher_data)} patients")
    axs[2].set_xticklabels(tmp_data["Gene"], rotation="vertical")
    axs[2].set_xlabel("")

    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
