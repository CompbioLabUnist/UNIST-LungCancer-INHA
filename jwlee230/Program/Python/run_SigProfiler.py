"""
run_SigProfiler.py: run SigProfilerExtractor as https://github.com/AlexandrovLab/SigProfilerExtractor
"""
import argparse
import os
import os.path
import shutil
import pandas
from SigProfilerExtractor import sigpro as sig
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input input Mutect2 VCF file(s)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data data CSV file", type=str)
    parser.add_argument("output", help="Output directory", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_histology = parser.add_mutually_exclusive_group(required=True)
    group_histology.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_histology.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".vcf"), args.input)):
        raise ValueError("Input files must end with .vcf!!")
    if not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be a positive integer!!")

    clinical_data: pandas.DataFrame = step00.get_clinical_data(args.clinical)
    clinical_types = sorted(clinical_data.columns)
    print(clinical_data)
    print(clinical_types)

    patients = list(map(step00.get_patient, args.input))

    if args.SQC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
        patients = list(filter(lambda x: step00.get_patient(x) in histology, patients))
    elif args.ADC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
        patients = list(filter(lambda x: step00.get_patient(x) in histology, patients))
    else:
        raise Exception("Something went wrong!!")
    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))
    print(sorted(patients))
    print(args.input)

    input_directory = os.path.join(args.output, "input")
    if os.path.exists(input_directory):
        shutil.rmtree(input_directory)

    os.makedirs(input_directory, mode=775)

    for input_file in args.input:
        os.system("cp -lv {0} {1}".format(input_file, os.path.join(input_directory, step00.get_id(input_file) + ".vcf")))

    sig.sigProfilerExtractor("vcf", args.output, input_data=args.output, reference_genome="GRCh38", opportunity_genome="GRCh38", exome=True, cpu=args.cpus)

    shutil.rmtree(input_directory)
