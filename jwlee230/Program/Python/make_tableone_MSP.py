"""
make_tableone.py: Create Table 1 for clinical data
"""
import argparse
import tableone
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".csv"):
        raise ValueError("Output must end with .CSV!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)
    print(sorted(clinical_data.columns))

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(patients)

    continuous_columns = ["Age", "Diagosis_Size_CT", "Diagosis_Size_OP", "Overall Survival", "Recurrence-Free Survival", "Pack-Year", "PSM"]
    categorical_columns = ["Gender", "Smoking-Detail", "pT", "pN", "Stage", "Vascular_invasion", "lymphatic_invasion", "Recurrence"]
    order = {"Gender": ["Male", "Female"], "Smoking-Detail": ["Never", "Ex", "Current"], "pT": [1, 2, 3, 4], "pN": [0, 1, 2, 3], "Stage": [1, 2, 3], "Vascular_invasion": ["NO", "YES"], "lymphatic_invasion": ["NO", "YES"], "Recurrence": ["NO", "YES"]}
    groupby = "PSM class"

    output_data = tableone.TableOne(clinical_data, columns=continuous_columns + categorical_columns, categorical=categorical_columns,  order=order, groupby=groupby, pval=True, htest_name=True, dip_test=True, normal_test=True, tukey_test=True, decimals={"PSM": 3, "PSM class": 0})
    print(output_data)
    output_data.to_csv(args.output)
