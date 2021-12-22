"""
make_input_mobster.py: make(select) input BAM file for mobster
"""
import argparse
import pandas
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input BAM files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data data CSV file", type=str)
    parser.add_argument("output", help="Output TXT file", type=str)

    group_histology = parser.add_mutually_exclusive_group(required=True)
    group_histology.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_histology.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".bam"), args.input)):
        raise ValueError("INPUT must end with .BAM!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".txt"):
        raise ValueError("Output must end with .TXT!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    patients = list(set(map(step00.get_patient, args.input)))
    if args.SQC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
        patients = list(filter(lambda x: step00.get_patient(x) in histology, patients))
    elif args.ADC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
        patients = list(filter(lambda x: step00.get_patient(x) in histology, patients))
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))
    args.input.sort(key=step00.sorting)

    with open(args.output, "w") as f:
        f.write(",".join(args.input))
