"""
calculate_CNV_DEG.py: calculate CNV for DEG
"""
import argparse
import logging
import typing
import gtfparse
import numpy
import pandas
import tqdm
import tqdm.contrib.concurrent

watching = ""
CNV_data = pandas.DataFrame()
gene_data = pandas.DataFrame()


def query(sample: str, chromosome: str, start: int, end: int) -> float:
    length = end - start + 1
    a = list()
    weights = list()

    tmp_data = CNV_data.loc[(CNV_data["Sample"] == sample) & (CNV_data["chromosome"] == chromosome) & (CNV_data["start"] <= start) & (end <= CNV_data["end"]), :]
    for index, row in tmp_data.iterrows():
        a.append(row[watching])
        weights.append(1)

    tmp_data = CNV_data.loc[(CNV_data["Sample"] == sample) & (CNV_data["chromosome"] == chromosome) & (start <= CNV_data["start"]) & (CNV_data["end"] <= end), :]
    for index, row in tmp_data.iterrows():
        a.append(row[watching])
        weights.append((row["end"] - row["start"] + 1) / length)

    tmp_data = CNV_data.loc[(CNV_data["Sample"] == sample) & (CNV_data["chromosome"] == chromosome) & (CNV_data["start"] <= start) & (start <= CNV_data["end"]), :]
    for index, row in tmp_data.iterrows():
        a.append(row[watching])
        weights.append((row["end"] - start + 1) / length)

    tmp_data = CNV_data.loc[(CNV_data["Sample"] == sample) & (CNV_data["chromosome"] == chromosome) & (CNV_data["start"] <= end) & (end <= CNV_data["end"]), :]
    for index, row in tmp_data.iterrows():
        a.append(row[watching])
        weights.append((end - row["start"] + 1) / length)

    if a and weights:
        return numpy.average(a=a, weights=weights)
    else:
        return 1.0


def run_gene(sample_and_gene: typing.Tuple[str, str]) -> float:
    sample, gene = sample_and_gene

    selected_gene_data = gene_data[(gene_data["gene_id"] == gene)]
    chromosome = list(selected_gene_data["seqname"])[0]
    start = min(selected_gene_data["start"])
    end = max(selected_gene_data["end"])

    return query(sample, chromosome, start, end)


if __name__ == "__main__":
    logging.getLogger("fontTools.subset").setLevel(logging.WARNING)

    parser = argparse.ArgumentParser()

    parser.add_argument("CNV", help="CNV segment.tsv file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("gene", help="Gene GTF file", type=str)
    parser.add_argument("cgc", help="CGC CSV file", type=str)
    parser.add_argument("output", help="Output TSV.gz file", type=str)
    parser.add_argument("--watching", help="Watching column name", type=str, required=True)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.CNV.endswith(".tsv"):
        raise ValueError("CNV must end with .TSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.gene.endswith(".gtf"):
        raise ValueError("Gene must end with .GTF!!")
    elif not args.cgc.endswith(".csv"):
        raise ValueError("CGC must end with .CSV!!")
    elif not args.output.endswith(".tsv.gz"):
        raise ValueError("Output must end with .TSV.gz!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    watching = args.watching

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(len(patients), patients)

    CNV_data = pandas.read_csv(args.CNV, sep="\t", index_col=0)
    CNV_data = CNV_data.loc[(CNV_data["Patient"].isin(patients))]
    print(CNV_data)

    sample_list = sorted(set(CNV_data["Sample"]))
    print("Sample:", len(sample_list))

    gene_data = gtfparse.read_gtf(args.gene)
    gene_data = gene_data.loc[(gene_data["feature"] == "exon")]
    print(gene_data)

    cgc_data = pandas.read_csv(args.cgc, index_col=0)
    print(cgc_data)

    gene_list = sorted(set(cgc_data.index) & set(gene_data["gene_id"]))
    # gene_list = sorted(set(gene_data["gene_id"]))
    print("Gene:", len(gene_list))

    output_data = pandas.DataFrame(index=gene_list, columns=sample_list, dtype=float)
    for sample in tqdm.tqdm(sample_list):
        output_data[sample] = tqdm.contrib.concurrent.process_map(run_gene, [(sample, gene) for gene in gene_list], max_workers=args.cpus, chunksize=1, leave=False)
    print(output_data)
    output_data.to_csv(args.output, sep="\t")
