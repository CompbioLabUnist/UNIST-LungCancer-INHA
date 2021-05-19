"""
make_clinical_MesKit.py: make clinical data for MesKit
"""
import argparse
import pandas
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input ID file(s)", type=str, nargs="+")
    parser.add_argument("suffix", help="Suffix", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)

    args = parser.parse_args()

    if not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")

    output_data = pandas.DataFrame()

    output_data["Tumor_Sample_Barcode"] = list(map(lambda x: x + args.suffix, args.input))
    output_data["Tumor_ID"] = list(map(step00.get_sample_type, output_data["Tumor_Sample_Barcode"]))
    output_data["Patient_ID"] = list(map(step00.get_patient, output_data["Tumor_Sample_Barcode"]))
    output_data["Tumor_Sample_Label"] = list(map(step00.get_long_sample_type, output_data["Tumor_Sample_Barcode"]))

    print(output_data)
    output_data.to_csv(args.output, sep="\t", index=False)
