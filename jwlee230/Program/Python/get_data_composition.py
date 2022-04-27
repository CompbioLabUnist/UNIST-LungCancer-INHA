"""
get_data_composition.py: get data composition by condition
"""
import argparse
import itertools
import pandas
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input IDs", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data data CSV file", type=str)
    parser.add_argument("output", help="Output file basename", type=str)

    args = parser.parse_args()

    if not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".tex"):
        raise ValueError("Output must end with .tex!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(sorted(clinical_data.columns))
    print(clinical_data)

    patients = list(map(step00.get_patient, args.input))
    print(clinical_data.loc[(clinical_data.index.isin(patients)) & (clinical_data["Histology"] == "SQC")])
    print(clinical_data.loc[(clinical_data.index.isin(patients)) & (clinical_data["Histology"] == "ADC")])

    output_data = pandas.DataFrame(columns=["Cancer Subtype", "Stage", "Number of Samples"], data=itertools.product(["LUSC", "LUAD"], step00.long_sample_type_list + ["Total"], [None])).set_index(keys=["Cancer Subtype", "Stage"], verify_integrity=True)
    print(output_data)

    for sample in tqdm.tqdm(args.input):
        raw_cancer_subtype = clinical_data.loc[step00.get_patient(sample), "Histology"]
        if raw_cancer_subtype == "SQC":
            cancer_subtype = "LUSC"
        elif raw_cancer_subtype == "ADC":
            cancer_subtype = "LUAD"
        else:
            raise ValueError("Unknown cancer subtype!!")

        stage = step00.get_long_sample_type(sample)

        if output_data.loc[(cancer_subtype, stage), "Number of Samples"] is None:
            output_data.loc[(cancer_subtype, stage), "Number of Samples"] = 0
        output_data.loc[(cancer_subtype, stage), "Number of Samples"] += 1

        if output_data.loc[(cancer_subtype, "Total"), "Number of Samples"] is None:
            output_data.loc[(cancer_subtype, "Total"), "Number of Samples"] = 0
        output_data.loc[(cancer_subtype, "Total"), "Number of Samples"] += 1

    output_data.dropna(inplace=True)
    print(output_data)

    output_data.to_latex(args.output, multirow=True, column_format="l|lr")
