"""
make_list_MutEnricher.py: Make VCF file list for MutEnricher input
"""
import argparse
import csv
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input VCF(.GZ) files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: (not x.endswith(".vcf")) and (not x.endswith(".vcf.gz")), args.input)):
        raise ValueError("INPUT must end with .VCF(.GZ)!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("OUTPUT must end with .TSV!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    args.input = list(filter(lambda x: step00.get_patient(step00.get_id(x)) in patients, args.input))
    args.input.sort()
    print(args.input)

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows([(i, step00.get_id(i)) for i in args.input])
