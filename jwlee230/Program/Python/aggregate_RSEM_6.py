"""
aggregate_RSEM_5.py: aggregate RSEM results for Mutation Shared Proportion
"""
import argparse
import datetime
import numpy
import pandas
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input RSEM filtered TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output file basename", type=str)
    parser.add_argument("--compare", help="Comparison MSP", type=str, choices=step00.sharing_columns)
    parser.add_argument("--date", help="Selection is date", action="store_true", default=False)

    group_histology = parser.add_mutually_exclusive_group(required=True)
    group_histology.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_histology.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_filter = parser.add_mutually_exclusive_group()
    group_filter.add_argument("--exact", help="Select only exact match (type, value)", type=str, nargs=2)
    group_filter.add_argument("--range", help="Select only in range (type, min, max)", type=str, nargs=3)

    group_threshold = parser.add_mutually_exclusive_group(required=True)
    group_threshold.add_argument("--median", help="Use median threshold", action="store_true", default=False)
    group_threshold.add_argument("--mean", help="Use mean threshold", action="store_true", default=False)

    args = parser.parse_args()

    if not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    clinical_types = sorted(clinical_data.columns)
    print(clinical_data)
    print(clinical_types)

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .TSV!!")
    elif (args.exact is not None) and ((args.exact[0] != "stage") and (args.exact[0] not in clinical_types)):
        raise ValueError("Exact should be one of the clinical types or <stage (sample type)>!!")
    elif (args.range is not None) and (args.range[0] not in clinical_types):
        raise ValueError("Range should be one of the clinical types!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col="gene_name")
    samples = sorted(input_data.columns, key=step00.sorting_by_type)
    print(input_data)
    print(samples)

    if args.SQC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
        samples = list(filter(lambda x: step00.get_patient(x) in histology, samples))
    elif args.ADC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
        samples = list(filter(lambda x: step00.get_patient(x) in histology, samples))
    else:
        raise Exception("Something went wrong!!")
    print(samples)

    if args.exact is not None:
        if args.exact[0] == "stage":
            assert args.exact[1] in step00.long_sample_type_list, "Invalid stage type!!"
            samples = list(filter(lambda x: step00.get_long_sample_type(x) == args.exact[1], samples))
        else:
            if args.date:
                value = datetime.datetime.strptime(args.exact[1], "%Y-%m-%d").date()
            else:
                value = args.exact[1]

            filtering = set(clinical_data.loc[(clinical_data[args.exact[0]] == value)].index)
            samples = list(filter(lambda x: step00.get_patient(x) in filtering, samples))
            assert samples, "Too less samples!!"
            print(samples)

    if args.range is not None:
        if args.date:
            min_value = datetime.datetime.strptime(args.range[1], "%Y-%m-%d").date()
            max_value = datetime.datetime.strptime(args.range[2], "%Y-%m-%d").date()
        else:
            min_value, max_value = args.range[1:]

        filtering = set(clinical_data.loc[(min_value <= clinical_data[args.range[0]]) & (clinical_data[args.range[0]] <= max_value)])
        samples = list(filter(lambda x: step00.get_patient(x) in filtering, samples))
        assert samples, "Too less samples!!"
        print(samples)

    if args.median:
        threshold = numpy.median(clinical_data.loc[list(set(map(step00.get_patient, samples))), args.compare])
    elif args.mean:
        threshold = numpy.mean(clinical_data.loc[list(set(map(step00.get_patient, samples))), args.compare])
    else:
        raise Exception("Something went wrong!!")
    print(threshold)

    clinical_data[f"{args.compare}-Lower/Higher"] = list(map(lambda x: "Lower" if (x <= threshold) else "Higher", clinical_data[args.compare]))
    print(clinical_data)

    f1 = set(clinical_data.loc[(clinical_data[f"{args.compare}-Lower/Higher"] == "Lower")].index)
    f2 = set(clinical_data.loc[(clinical_data[f"{args.compare}-Lower/Higher"] == "Higher")].index)

    samples = list(filter(lambda x: step00.get_patient(x) in f1, samples)) + list(filter(lambda x: step00.get_patient(x) in f2, samples))
    conditions = list(map(lambda x: "Lower" if (step00.get_patient(x) in f1) else "Higher", samples))

    assert samples, "Too less samples!!"
    assert len(set(conditions)) == 2, "Invalid conditions!!"

    print(samples)
    print(conditions)

    input_data.loc[:, samples].to_csv(args.output + ".tsv", sep="\t", index=True, header=True)
    pandas.DataFrame(data=zip(samples, conditions), columns=["ID", "condition"]).to_csv(args.output + ".coldata", sep="\t", index=False, header=True)
