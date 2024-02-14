"""
calculate_DEG_MSP-ratio.py: calculate correlation between DEG and MSP-ratio
"""
import argparse
import multiprocessing
import numpy
import pandas
import scipy.stats
import tqdm
import step00

input_data = pandas.DataFrame()
clinical_data = pandas.DataFrame()


def correlation(MSP: str, gene: str):
    precancer_list = list(filter(lambda x: x in set(input_data.index), clinical_data[f"{MSP}-sample"]))
    primary_list = list(map(step00.get_paired_primary, precancer_list))

    MSP_list = input_data.loc[precancer_list, MSP]
    precancer = numpy.array(input_data.loc[precancer_list, gene]) + step00.epsilon
    primary = numpy.array(input_data.loc[primary_list, gene]) + step00.epsilon
    ratio_list = numpy.divide(precancer, primary)

    r = scipy.stats.linregress(MSP_list, ratio_list)
    return list(r) + [r.intercept_stderr]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input gene TSV file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with.TSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0).T
    genes = list(input_data.columns)
    print(input_data)

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

    input_data = input_data.loc[list(filter(lambda x: step00.get_patient(x) in patients, list(input_data.index))), :]
    input_data["Stage"] = list(map(step00.get_long_sample_type, list(input_data.index)))
    for column in tqdm.tqdm(step00.sharing_columns):
        input_data[column] = list(map(lambda x: clinical_data.loc[step00.get_patient(x), column], list(input_data.index)))
    print(input_data)

    columns = ["slope", "intercept", "r", "p", "stderr", "intercept_stderr"]

    output_data = pandas.DataFrame(index=genes)
    with multiprocessing.Pool(args.cpus) as pool:
        for MSP in tqdm.tqdm(step00.sharing_columns):
            output_data = output_data.join(pandas.DataFrame(pool.starmap(correlation, [(MSP, gene) for gene in genes]), index=genes, columns=[f"{MSP}-{column}" for column in columns]))

            output_data[f"{MSP}-log10(abs(slope))"] = list(map(lambda x: numpy.log10(abs(x)) if x != 0 else 0, output_data[f"{MSP}-slope"]))
            output_data[f"{MSP}-importance"] = list(map(lambda x: abs(x[0] * x[1]), zip(output_data[f"{MSP}-r"], output_data[f"{MSP}-slope"])))

    print(output_data)
    output_data.to_csv(args.output, sep="\t")
