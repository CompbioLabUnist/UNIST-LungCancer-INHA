"""
merge_CNVkit_Gistic_2.py: merge CNVkit CNS files for Gistic with cancer stage
"""
import argparse
import multiprocessing
import pandas
import tqdm
import step00


def get_data(file_name: str) -> pandas.DataFrame:
    data = pandas.read_csv(file_name, sep="\t")
    data["ID"] = step00.get_id(file_name)
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PureCN output CNS file(s)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("blacklist", help="UCSC blacklist file", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--stage", help="Stage choices", choices=step00.long_sample_type_list, nargs="+", required=True)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".cns"), args.input)):
        raise ValueError("INPUT must end with .CNS!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    blacklist_data = pandas.read_csv(args.blacklist, sep="\t")
    print(blacklist_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(patients)
    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))

    print(args.stage)
    args.input = list(filter(lambda x: step00.get_long_sample_type(x) in args.stage, args.input))

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="index", copy=False, ignore_index=True, verify_integrity=True)
    print(input_data)

    for index, row in tqdm.tqdm(blacklist_data.iterrows()):
        chrom, start, end = row["chrom"], row["chromStart"], row["chromEnd"]
        input_data = input_data.loc[~((input_data["chromosome"] == chrom) & (start <= input_data["start"]) & (input_data["end"] <= end)), :]

    input_data = input_data.loc[:, ["ID", "chromosome", "start", "end", "probes", "log2"]]
    input_data.columns = step00.Gistic_columns
    print(input_data)

    input_data.to_csv(args.output, sep="\t", index=False)
