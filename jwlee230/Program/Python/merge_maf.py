"""
merge_maf.py: merge MAF files for sharing
"""
import argparse
import multiprocessing
import pandas
import step00


def read_file(filename):
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 MAF file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("Input must end with .MAF!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    input_data = pandas.DataFrame(index=list(map(lambda x: step00.get_id(x), args.input)))
    with multiprocessing.Pool(args.cpus) as pool:
        output_data = pandas.concat(pool.map(read_file, args.input), ignore_index=True)
        output_data = output_data.loc[(output_data[step00.nonsynonymous_column].isin(step00.nonsynonymous_mutations))]
        output_data["Tumor_Sample_Barcode"] = list(map(step00.get_id, output_data["Tumor_Sample_Barcode"]))
        output_data["Matched_Norm_Sample_Barcode"] = list(map(step00.get_id, output_data["Matched_Norm_Sample_Barcode"]))

    output_data = output_data[["Hugo_Symbol", "Entrez_Gene_Id", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2"]]
    print(output_data)
    output_data.to_csv(args.output, sep="\t", index=False)
