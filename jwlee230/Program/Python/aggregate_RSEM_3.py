"""
aggregate_RSEM_3.py: aggregate RSEM results for recurrence/non-recurrence with merging normal samples
"""
import argparse
import datetime
import pandas
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input RSEM filtered TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data data CSV file", type=str)
    parser.add_argument("output", help="Output file basename", type=str)
    parser.add_argument("--compare", help="Comparison grouping (type, control, case)", type=str, nargs=3, default=["stage", "Normal", "Primary"])
    parser.add_argument("--date", help="Selection is date", action="store_true", default=False)

    group_histology = parser.add_mutually_exclusive_group(required=True)
    group_histology.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_histology.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_filter = parser.add_mutually_exclusive_group()
    group_filter.add_argument("--exact", help="Select only exact match (type, value)", type=str, nargs=2)
    group_filter.add_argument("--range", help="Select only in range (type, min, max)", type=str, nargs=3)

    args = parser.parse_args()

    if not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    clinical_types = sorted(clinical_data.columns)
    print(clinical_data)
    print(clinical_types)

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .TSV!!")
    elif (args.compare[0] != "stage") and (args.compare[0] not in clinical_types):
        raise ValueError("Campare should be one of the clinical types or <stage (sample type)>!!")
    elif (args.exact is not None) and ((args.exact[0] != "stage") and (args.exact[0] not in clinical_types)):
        raise ValueError("Exact should be one of the clinical types or <stage (sample type)>!!")
    elif (args.range is not None) and (args.range[0] not in clinical_types):
        raise ValueError("Range should be one of the clinical types!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col="gene_name")
    original_patients = sorted(input_data.columns, key=step00.sorting_by_type)
    print(input_data)
    print(original_patients)

    if args.SQC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
        original_patients = list(filter(lambda x: step00.get_patient(x) in histology, original_patients))
    elif args.ADC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
        original_patients = list(filter(lambda x: step00.get_patient(x) in histology, original_patients))
    else:
        raise Exception("Something went wrong!!")
    print(original_patients)

    if args.exact is not None:
        if args.exact[0] == "stage":
            assert args.exact[1] in step00.long_sample_type_list, "Invalid stage type!!"
            patients = list(filter(lambda x: step00.get_long_sample_type(x) == args.exact[1], original_patients))
        else:
            if args.date:
                value = datetime.datetime.strptime(args.exact[1], "%Y-%m-%d").date()
            else:
                value = args.exact[1]

            filtering = set(clinical_data.loc[(clinical_data[args.exact[0]] == value)].index)
            patients = list(filter(lambda x: step00.get_patient(x) in filtering, original_patients))
            assert patients, "Too less patients!!"
            print(patients)

    if args.range is not None:
        if args.date:
            min_value = datetime.datetime.strptime(args.range[1], "%Y-%m-%d").date()
            max_value = datetime.datetime.strptime(args.range[2], "%Y-%m-%d").date()
        else:
            min_value, max_value = args.range[1:]

        filtering = set(clinical_data.loc[(min_value <= clinical_data[args.range[0]]) & (clinical_data[args.range[0]] <= max_value)])
        patients = list(filter(lambda x: step00.get_patient(x) in filtering, patients))
        assert patients, "Too less patients!!"
        print(patients)

    if args.compare[0] == "stage":
        assert args.compare[1] in step00.long_sample_type_list, "Invalid stage type!!"
        assert args.compare[2] in step00.long_sample_type_list, "Invalid stage type!!"

        patients = list(filter(lambda x: step00.get_long_sample_type(x) == args.compare[1], original_patients if (args.compare[1] == "Normal") else patients)) + list(filter(lambda x: step00.get_long_sample_type(x) == args.compare[2], original_patients if (args.compare[2] == "Normal") else patients))
        conditions = list(map(step00.get_long_sample_type, patients))

        assert patients, "Too less patients!!"
        assert len(set(conditions)) == 2, "Invalid conditions!!"

        print(patients)
        print(conditions)
    else:
        f1 = set(clinical_data.loc[(clinical_data[args.compare[0]] == args.compare[1])].index)
        f2 = set(clinical_data.loc[(clinical_data[args.compare[0]] == args.compare[2])].index)

        patients = list(filter(lambda x: step00.get_patient(x) in f1, patients)) + list(filter(lambda x: step00.get_patient(x) in f2, patients))
        conditions = list(map(lambda x: args.compare[1] if (step00.get_patient(x) in f1) else args.compare[2], patients))

        assert patients, "Too less patients!!"
        assert len(set(conditions)) == 2, "Invalid conditions!!"

        print(patients)
        print(conditions)

    input_data.loc[:, patients].to_csv(args.output + ".tsv", sep="\t", index=True, header=True)
    pandas.DataFrame(data=zip(patients, conditions), columns=["ID", "condition"]).to_csv(args.output + ".coldata", sep="\t", index=False, header=True)
