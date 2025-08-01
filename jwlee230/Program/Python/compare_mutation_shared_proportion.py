"""
compare_mutation_shared_proportion.py: Compare mutation shared proportion
"""
import argparse
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import seaborn
import pandas
import tqdm
import step00

mutect_data = pandas.DataFrame()


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


def query_mutation(gene: str, patient: str) -> int:
    patient_data = mutect_data.loc[(mutect_data["Patient"] == patient) & (mutect_data["Hugo_Symbol"] == gene)]

    stage_set = list(filter(lambda x: x in set(patient_data["Stage"]), step00.long_sample_type_list))
    if ("Primary" in stage_set) and (len(stage_set) > 1):
        primary_set = set(patient_data.loc[patient_data["Stage"] == "Primary", step00.sharing_strategy].itertuples(index=False, name=None))
        precancer_set = set(patient_data.loc[patient_data["Stage"] == stage_set[-2], step00.sharing_strategy].itertuples(index=False, name=None))
        return len(precancer_set & primary_set)
    else:
        return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion", type=str)
    parser.add_argument("cgc", help="CGC gene CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--threshold", help="Threshold to use", type=int, default=30)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_strategy = parser.add_mutually_exclusive_group(required=True)
    group_strategy.add_argument("--median", help="Median division", action="store_true", default=False)
    group_strategy.add_argument("--mean", help="Mean division", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.cgc.endswith(".csv"):
        raise ValueError("CGC must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif args.threshold <= 5:
        raise ValueError("Threshold must be greater than 5!!")

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(patients)

    cgc_data = pandas.read_csv(args.cgc, index_col="Gene Symbol")
    print(cgc_data)

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))
    args.input.sort(key=step00.sorting)
    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
        mutect_data["Tumor_Sample_Barcode"] = pool.map(step00.get_id, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Patient"] = pool.map(step00.get_patient, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Stage"] = pool.map(step00.get_long_sample_type, mutect_data["Tumor_Sample_Barcode"])
    print(mutect_data)

    mutect_data = mutect_data.loc[(mutect_data[step00.nonsynonymous_column].isin(step00.nonsynonymous_mutations))]
    print(mutect_data)

    patients &= set(mutect_data["Patient"])

    figures = list()
    for column in step00.sharing_columns:
        if args.median:
            cutting = numpy.median(clinical_data[column])
        elif args.mean:
            cutting = numpy.mean(clinical_data[column])
        else:
            raise Exception("Something went wrong!!")

        lower_data = clinical_data.loc[(clinical_data[column] <= cutting)]
        higher_data = clinical_data.loc[(clinical_data[column] > cutting)]
        print(lower_data)
        print(higher_data)

        data = list()
        with multiprocessing.Pool(args.cpus) as pool:
            for gene in tqdm.tqdm(list(cgc_data.index)):
                data.append((gene, sum(pool.starmap(query_mutation, [(gene, patient) for patient in list(lower_data.index)])), "Lower"))
                data.append((gene, sum(pool.starmap(query_mutation, [(gene, patient) for patient in list(higher_data.index)])), "Higher"))
        output_data = pandas.DataFrame(data=data, columns=["Gene", "Count", "Lower/Higher"])
        print(output_data)

        matplotlib.use("Agg")
        matplotlib.rcParams.update(step00.matplotlib_parameters)
        seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

        fig, axs = matplotlib.pyplot.subplots(figsize=(32 * 3, 18), ncols=3, sharey=True)

        genes = sorted(list(cgc_data.index), key=lambda x: sum(output_data.loc[(output_data["Gene"] == x), "Count"]), reverse=True)[:args.threshold]
        tmp_data = output_data.loc[(output_data["Gene"].isin(genes))].sort_values(by="Count", ascending=False)
        seaborn.barplot(data=tmp_data, x="Gene", y="Count", order=genes, hue="Lower/Higher", ci=None, palette={"Lower": "tab:cyan", "Higher": "tab:red"}, ax=axs[0])
        axs[0].set_title(f"Total {len(clinical_data)} patients")
        axs[0].set_xticklabels(genes, rotation="vertical")
        axs[0].set_xlabel("")
        axs[0].set_ylabel("Patients")

        tmp_data = output_data.loc[(output_data["Lower/Higher"] == "Lower")].sort_values(by="Count", ascending=False).iloc[:args.threshold, :]
        seaborn.barplot(data=tmp_data, x="Gene", y="Count", color="tab:cyan", ci=None, ax=axs[1])
        axs[1].set_title(f"Lower {len(lower_data)} patients")
        axs[1].set_xticklabels(tmp_data["Gene"], rotation="vertical")
        axs[1].set_xlabel("")
        axs[1].set_ylabel("Patients")

        tmp_data = output_data.loc[(output_data["Lower/Higher"] == "Higher")].sort_values(by="Count", ascending=False).iloc[:args.threshold, :]
        seaborn.barplot(data=tmp_data, x="Gene", y="Count", color="tab:red", ci=None, ax=axs[2])
        axs[2].set_title(f"Higher {len(higher_data)} patients")
        axs[2].set_xticklabels(tmp_data["Gene"], rotation="vertical")
        axs[2].set_xlabel("")
        axs[2].set_ylabel("Patients")

        matplotlib.pyplot.tight_layout()

        figures.append(f"{column}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
