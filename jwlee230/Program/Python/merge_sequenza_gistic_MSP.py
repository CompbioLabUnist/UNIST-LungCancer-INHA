"""
merge_sequenza_gistic_MSP.py: merge Sequenza segment file for GISTIC
"""
import argparse
import multiprocessing
import numpy
import pandas
import tqdm
import step00


def get_data(file_name: str) -> pandas.DataFrame:
    data = pandas.read_csv(file_name, sep="\t")
    data["ID"] = step00.get_id(file_name.split("/")[-2])
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Sequenza output segment TXT file(s)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data w/ MSP TSV file", type=str)
    parser.add_argument("blacklist", help="UCSC blacklist file", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)
    parser.add_argument("--threshold", help="Threshold for MSP", type=int, default=25)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_updown = parser.add_mutually_exclusive_group(required=True)
    group_updown.add_argument("--up", help="Get segment with MSP >= threshold", action="store_true", default=False)
    group_updown.add_argument("--down", help="Get segment with MSP <= threshold", action="store_true", default=False)

    group_precancer = parser.add_mutually_exclusive_group(required=True)
    group_precancer.add_argument("--pre", help="Get precancer sample only", action="store_true", default=False)
    group_precancer.add_argument("--pri", help="Get primary tumor sample only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".txt"), args.input)):
        raise ValueError("INPUT must end with .txt!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    print(clinical_data)

    threshold = numpy.percentile(clinical_data[step00.sharing_columns[1]], args.threshold)
    if args.up:
        clinical_data = clinical_data.loc[(clinical_data[step00.sharing_columns[1]] >= threshold)]
    elif args.down:
        clinical_data = clinical_data.loc[(clinical_data[step00.sharing_columns[1]] <= threshold)]
    else:
        raise Exception("Something went wrong!!")
    print(clinical_data)

    blacklist_data = pandas.read_csv(args.blacklist, sep="\t")
    print(blacklist_data)

    if args.pre:
        args.input = list(filter(lambda x: x.split("/")[-2] in set(clinical_data[step00.sharing_columns[1] + "-sample"]), args.input))
    elif args.pri:
        args.input = list(filter(lambda x: step00.get_id(x.split("/")[-2]) in set(map(step00.get_paired_primary, clinical_data[step00.sharing_columns[1] + "-sample"])), args.input))
    else:
        raise Exception("Something went wrong!!")

    sample_list = list(map(lambda x: step00.get_id(x.split("/")[-2]), args.input))
    print(len(sample_list), sample_list)

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="index", copy=False, ignore_index=True, verify_integrity=True)
    print(input_data)

    for index, row in tqdm.tqdm(blacklist_data.iterrows(), total=len(blacklist_data)):
        chrom, start, end = row["chrom"], row["chromStart"], row["chromEnd"]
        input_data = input_data.loc[~((input_data["chromosome"] == chrom) & (start <= input_data["start.pos"]) & (input_data["end.pos"] <= end)), :]

    input_data = input_data.loc[:, ["ID", "chromosome", "start.pos", "end.pos", "N.BAF", "depth.ratio"]]
    input_data.columns = step00.Gistic_columns
    print(input_data)
    input_data.to_csv(args.output, sep="\t", index=False)
