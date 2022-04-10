"""
aggregate_mutect_deconstructSigs.py: Aggregate Mutect2 MAF results with deconstructSigs
"""
import argparse
import multiprocessing
import pandas
import step00

columns = ["Tumor_Sample_Barcode", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"]


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False, usecols=columns)[columns]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .tsv!!")
    elif args.cpus < 1:
        raise ValueError("CPUS must be positive!!")

    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
        mutect_data = mutect_data.loc[(mutect_data["Chromosome"].isin(step00.chromosome_list)) & (mutect_data["Reference_Allele"].str.len() == 1) & (mutect_data["Tumor_Seq_Allele2"].str.len() == 1), :]
        mutect_data["Tumor_Sample_Barcode"] = pool.map(step00.get_id, mutect_data["Tumor_Sample_Barcode"])
    mutect_data.columns = ["Sample", "chr", "pos", "ref", "alt"]
    print(mutect_data)

    mutect_data.to_csv(args.output, sep="\t", index=False)
