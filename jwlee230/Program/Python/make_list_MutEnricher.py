"""
make_list_MutEnricher.py: Make VCF file list for MutEnricher input
"""
import argparse
import csv

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input VCF(.GZ) files", type=str, nargs="+")
    parser.add_argument("output", help="Output TSV file", type=str)

    args = parser.parse_args()

    if list(filter(lambda x: (not x.endswith(".vcf")) and (not x.endswith(".vcf.gz")), args.input)):
        raise ValueError("INPUT must end with .VCF(.GZ)!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("OUTPUT must end with .TSV!!")

    args.input.sort()

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows([(i, i.split("/")[-1].split(".")[0]) for i in args.input])
