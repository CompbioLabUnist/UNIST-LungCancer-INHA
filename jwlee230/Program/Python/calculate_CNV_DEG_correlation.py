"""
calculate_CNV_DEG_correlation.py: calculate CNV & DEG correlation
"""
import argparse
import itertools
import multiprocessing
import typing
import numpy
import pandas
import scipy.stats
import tqdm
import step00

input_data = pandas.DataFrame()


def correlation_stage(stage: str, sharing: str, gene: str) -> typing.Tuple[float, float, float, float, float]:
    tmp_data = input_data.loc[(input_data["Stage"] == stage), [sharing, gene]]
    if tmp_data.shape[0] < 3:
        return 0.0, 0.0, 0.0, 1.0, 0.0
    elif (numpy.std(tmp_data[sharing]) == 0.0) or (numpy.std(tmp_data[gene]) == 0.0):
        return 0.0, 0.0, 0.0, 1.0, 0.0
    return tuple(scipy.stats.linregress(tmp_data[sharing], tmp_data[gene]))


def correlation_precancer(sharing: str, gene: str) -> typing.Tuple[float, float, float, float, float]:
    tmp_data = input_data.loc[~(input_data["Stage"].isin({"Normal", "Primary"})), [sharing, gene]]
    if tmp_data.shape[0] < 3:
        return 0.0, 0.0, 0.0, 1.0, 0.0
    elif (numpy.std(tmp_data[sharing]) == 0.0) or (numpy.std(tmp_data[gene]) == 0.0):
        return 0.0, 0.0, 0.0, 1.0, 0.0
    return tuple(scipy.stats.linregress(tmp_data[sharing], tmp_data[gene]))


def correlation_all(sharing: str, gene: str) -> typing.Tuple[float, float, float, float, float]:
    tmp_data = input_data.loc[:, [sharing, gene]]
    if tmp_data.shape[0] < 3:
        return 0.0, 0.0, 0.0, 1.0, 0.0
    elif (numpy.std(tmp_data[sharing]) == 0.0) or (numpy.std(tmp_data[gene]) == 0.0):
        return 0.0, 0.0, 0.0, 1.0, 0.0
    return tuple(scipy.stats.linregress(tmp_data[sharing], tmp_data[gene]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input gene TSV.gz file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv.gz"):
        raise ValueError("Input must end with .TSV.gz!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0).T
    genes = list(input_data.columns)
    print(input_data)

    clinical_data: pandas.DataFrame = pandas.read_csv(args.clinical, sep="\t", index_col=0)
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

    stages = list(filter(lambda x: x in set(input_data["Stage"]), step00.long_sample_type_list)) + ["Precancer", "All"]
    columns = ["slope", "intercept", "r", "p", "stderr"]

    output_data = pandas.DataFrame(index=genes)
    with multiprocessing.Pool(args.cpus) as pool:
        for MSP, stage in tqdm.tqdm(list(itertools.product(step00.sharing_columns, stages))):
            if stage == "All":
                output_data = output_data.join(pandas.DataFrame(pool.starmap(correlation_all, [(MSP, gene) for gene in genes]), index=genes, columns=[f"{stage}-{MSP}-{column}" for column in columns]))
            elif stage == "Precancer":
                output_data = output_data.join(pandas.DataFrame(pool.starmap(correlation_precancer, [(MSP, gene) for gene in genes]), index=genes, columns=[f"{stage}-{MSP}-{column}" for column in columns]))
            else:
                output_data = output_data.join(pandas.DataFrame(pool.starmap(correlation_stage, [(stage, MSP, gene) for gene in genes]), index=genes, columns=[f"{stage}-{MSP}-{column}" for column in columns]))

            output_data[f"{stage}-{MSP}-slope"] = list(map(abs, output_data[f"{stage}-{MSP}-slope"]))
            output_data[f"{stage}-{MSP}-log10(abs(slope))"] = list(map(lambda x: numpy.log10(x) if (x > step00.epsilon) else 0, output_data[f"{stage}-{MSP}-slope"]))
            output_data[f"{stage}-{MSP}-importance"] = list(map(lambda x: abs(x[0] * x[1]), zip(output_data[f"{stage}-{MSP}-r"], output_data[f"{stage}-{MSP}-slope"])))
    print(output_data)
    output_data.to_csv(args.output, sep="\t")
