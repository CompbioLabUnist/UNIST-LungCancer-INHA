"""
get_data_composition_gene.py: get data composition which separated with gene mutation
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
    parser.add_argument("--mutect", help="Mutect2 results MAF file(s)", type=str, nargs="+", required=True)
    parser.add_argument("--gene", help="Gene name for separate data", type=str, required=True)

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

    args.input.sort(key=step00.sorting_by_type)
    args.mutect.sort(key=step00.sorting_by_type)

    if args.SQC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
        args.input = list(filter(lambda x: step00.get_patient(x) in histology, args.input))
    elif args.ADC:
        histology = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
        args.input = list(filter(lambda x: step00.get_patient(x) in histology, args.input))
    else:
        raise Exception("Something went wrong!!")
    print(args.input)

    output_data = pandas.DataFrame(columns=[args.gene, "Stage", "Number of Samples"], data=itertools.product(["No", "Synonymous", "Non-synonymous"], step00.long_sample_type_list + ["Total"], [None]))
    output_data.loc[-1] = ["Normal", "Normal", 0]
    output_data = output_data.sort_index().set_index(keys=[args.gene, "Stage"], verify_integrity=True)
    print(output_data)

    for sample in tqdm.tqdm(args.input):
        stage = step00.get_long_sample_type(sample)

        if stage == "Normal":
            output_data.loc[("Normal", stage), "Number of Samples"] += 1
            continue

        maf_file = list(filter(lambda x: step00.get_id(x) == sample, args.mutect))[0]
        maf_data = pandas.read_csv(maf_file, sep="\t", comment="#", usecols=["Hugo_Symbol", "Variant_Classification"])

        maf_data = maf_data.loc[(maf_data["Hugo_Symbol"] == args.gene)]
        if maf_data.empty:
            if output_data.loc[("No", stage), "Number of Samples"] is None:
                output_data.loc[("No", stage), "Number of Samples"] = 0
            output_data.loc[("No", stage), "Number of Samples"] += 1
            continue

        maf_data = maf_data.loc[(maf_data["Variant_Classification"].isin(step00.nonsynonymous_mutations))]
        if maf_data.empty:
            if output_data.loc[("Synonymous", stage), "Number of Samples"] is None:
                output_data.loc[("Synonymous", stage), "Number of Samples"] = 0
            output_data.loc[("Synonymous", stage), "Number of Samples"] += 1
        else:
            if output_data.loc[("Non-synonymous", stage), "Number of Samples"] is None:
                output_data.loc[("Non-synonymous", stage), "Number of Samples"] = 0
            output_data.loc[("Non-synonymous", stage), "Number of Samples"] += 1

    output_data.dropna(inplace=True)
    print(output_data)

    output_data.to_latex(args.output, multirow=True, column_format="l|lr")
