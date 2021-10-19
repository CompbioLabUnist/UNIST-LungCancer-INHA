"""
get_data_composition_smoking.py: get data composition which separated with smoking condition
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

    group_histology = parser.add_mutually_exclusive_group(required=True)
    group_histology.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_histology.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".tex"):
        raise ValueError("Output must end with .tex!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(sorted(clinical_data.columns))
    print(clinical_data)

    if args.SQC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
        args.input = list(filter(lambda x: step00.get_patient(x) in histology, args.input))
    elif args.ADC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
        args.input = list(filter(lambda x: step00.get_patient(x) in histology, args.input))
    else:
        raise Exception("Something went wrong!!")
    print(args.input)

    output_data = pandas.DataFrame(columns=["Smoking?", "Stage", "Number of Samples"], data=itertools.product(["Never", "Ever"], step00.long_sample_type_list + ["Total"], [None])).set_index(keys=["Smoking?", "Stage"], verify_integrity=True)
    print(output_data)

    for sample in tqdm.tqdm(args.input):
        smoking = clinical_data.loc[step00.get_patient(sample), "Smoking"]
        stage = step00.get_long_sample_type(sample)

        if output_data.loc[(smoking, stage), "Number of Samples"] is None:
            output_data.loc[(smoking, stage), "Number of Samples"] = 0
        output_data.loc[(smoking, stage), "Number of Samples"] += 1

        if output_data.loc[(smoking, "Total"), "Number of Samples"] is None:
            output_data.loc[(smoking, "Total"), "Number of Samples"] = 0
        output_data.loc[(smoking, "Total"), "Number of Samples"] += 1

    output_data.dropna(inplace=True)
    print(output_data)

    output_data.to_latex(args.output, multirow=True, column_format="l|lr")
