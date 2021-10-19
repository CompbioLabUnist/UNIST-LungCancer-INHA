"""
make_sample_map_haplotypecaller.py: make sample map file for HaplotypeCaller
"""
import argparse
import csv
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input Haplotype Caller output g.vcf.gz file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--patient", help="Designated patient name", required=True, type=str)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".g.vcf.gz"), args.input)):
        raise ValueError("Input file must end in .g.vcf.gz!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output file must end in .tsv!!")

    args.input.sort(key=step00.sorting)
    sample_list = list(map(lambda x: (x.split("/")[-1].split(".")[0], x), list(filter(lambda x: step00.get_patient(x) == args.patient, args.input))))
    assert sample_list, "No patient remaining!!"

    with open(args.output, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(sample_list)
