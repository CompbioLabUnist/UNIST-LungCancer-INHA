"""
plot_gene_clinical_venn.py: Plot venn diagram the importance of gene upon clinical data
"""
import argparse
import multiprocessing
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import tqdm
import venn
import step00

mutect_data = pandas.DataFrame()


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


def query(gene: str, stage: str) -> int:
    return mutect_data.loc[(mutect_data["Cancer_subtype"] == stage) & (mutect_data["Hugo_Symbol"] == gene), :].shape[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("figure", help="Output PDF file", type=str)
    parser.add_argument("table", help="Output TSV files", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.figure.endswith(".pdf"):
        raise ValueError("Figure must end with .PDF!!")
    elif not args.table.endswith(".tsv"):
        raise ValueError("Table must end with .TSV!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(sorted(patients))

    args.input = sorted(filter(lambda x: step00.get_patient(x.split("/")[-1].split(".")[0]) in patients, args.input), key=step00.sorting_by_type)

    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)

    mutect_data = mutect_data.loc[(mutect_data["Variant_Classification"].isin(step00.nonsynonymous_mutations))]
    mutect_data["Tumor_Sample_Barcode"] = list(map(lambda x: x.split(".")[0], mutect_data["Tumor_Sample_Barcode"]))
    mutect_data["Cancer_subtype"] = list(map(step00.get_long_sample_type, mutect_data["Tumor_Sample_Barcode"]))
    print(mutect_data)

    selected_stage_list = list(filter(lambda x: x in set(mutect_data["Cancer_subtype"]), step00.long_sample_type_list))
    print(selected_stage_list)

    input_data = dict()
    for stage in tqdm.tqdm(selected_stage_list):
        input_data[stage] = set(mutect_data.loc[(mutect_data["Cancer_subtype"] == stage), "Hugo_Symbol"])

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))
    venn.venn(input_data, ax=ax, fmt=step00.venn_format, fontsize=step00.matplotlib_parameters["legend.fontsize"], legend_loc="upper left")

    matplotlib.pyplot.tight_layout()
    fig.savefig(args.figure)
    matplotlib.pyplot.close(fig)

    gene_set = set()
    for stage in tqdm.tqdm(selected_stage_list):
        gene_set |= input_data[stage]

    table_data = pandas.DataFrame(index=sorted(gene_set), columns=selected_stage_list, dtype=int)
    with multiprocessing.Pool(args.cpus) as pool:
        for stage in tqdm.tqdm(selected_stage_list):
            table_data.loc[:, stage] = pool.starmap(query, [(gene, stage) for gene in sorted(gene_set)])
    print(table_data)
    table_data.to_csv(args.table, sep="\t")
