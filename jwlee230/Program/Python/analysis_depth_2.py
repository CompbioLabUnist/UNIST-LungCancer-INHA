"""
Python/analysis_depth_2.py: analysis WES depth
"""
import argparse
import typing
import multiprocessing
import numpy
import pandas
import tqdm
import step00


def get_depth(file_name: str) -> typing.Tuple[str, float, float, int]:
    """
    get_depth: get depth of coverage from file
    """
    data = pandas.read_csv(file_name, sep="\t", names=["CHROM", "POS", "Depth"], skiprows=1)
    return (step00.get_id(file_name), numpy.mean(data["Depth"]), numpy.std(data["Depth"]) / 2, sum(data["Depth"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Samtools depth TSV files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data TSV file", type=str)
    parser.add_argument("output", help="Output TST file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".txt"):
        raise ValueError("Output must end with .TXT!!")
    elif args.cpus < 0:
        raise ValueError("CPUS must be positive!!")

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(patients)

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))

    with multiprocessing.Pool(processes=args.cpus) as pool:
        input_data = pandas.DataFrame(data=pool.map(get_depth, args.input), columns=["ID", "mean", "std", "sum"])
    print(input_data)

    for MSP in tqdm.tqdm(step00.sharing_columns):
        print(MSP)

        precancer_list = list(clinical_data[f"{MSP}-sample"])
        normal_list = list(map(step00.get_paired_normal, precancer_list))
        primary_list = list(map(step00.get_paired_primary, precancer_list))

        total_data = input_data.loc[(input_data["ID"].isin(precancer_list + normal_list + primary_list))]
        print("Reads:", numpy.mean(total_data["mean"]), numpy.std(total_data["mean"]))

        normal_data = input_data.loc[(input_data["ID"].isin(normal_list))]
        print("Normal:", numpy.mean(normal_data["mean"]), numpy.std(normal_data["mean"]))

        precancer_data = input_data.loc[(input_data["ID"].isin(precancer_list))]
        print("Precancer:", numpy.mean(precancer_data["mean"]), numpy.std(precancer_data["mean"]))

        primary_data = input_data.loc[(input_data["ID"].isin(primary_list))]
        print("Primary:", numpy.mean(primary_data["mean"]), numpy.std(primary_data["mean"]))

        print()
