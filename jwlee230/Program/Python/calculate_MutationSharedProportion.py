"""
calculate_MutationSharedProportion.py: calculate Mutation Shared Proportion
"""
import argparse
import multiprocessing
import numpy
import pandas
import tqdm
import step00


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF file(s)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    clinical_data: pandas.DataFrame = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    patients = set(clinical_data.index)
    print(len(patients))

    args.input.sort(key=step00.sorting)
    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
        mutect_data["Tumor_Sample_Barcode"] = pool.map(step00.get_id, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Patient"] = pool.map(step00.get_patient, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Stage"] = pool.map(step00.get_long_sample_type, mutect_data["Tumor_Sample_Barcode"])
    print(mutect_data)

    for MSP in tqdm.tqdm(step00.sharing_columns):
        clinical_data[MSP] = None
        clinical_data[f"{MSP}-sample"] = None

    for patient in tqdm.tqdm(patients):
        patient_data = mutect_data.loc[(mutect_data["Patient"] == patient) & (mutect_data[step00.nonsynonymous_column].isin(step00.nonsynonymous_mutations))]
        primary_set = set(patient_data.loc[(patient_data["Stage"] == "Primary"), step00.sharing_strategy].itertuples(index=False, name=None))
        precancer_list = list(mutect_data.loc[(mutect_data["Patient"] == patient) & ~(mutect_data["Stage"].isin({"Primary"})), "Tumor_Sample_Barcode"])

        if (not precancer_list) or (not primary_set):
            continue

        for precancer in precancer_list:
            precancer_set = set(patient_data.loc[(patient_data["Tumor_Sample_Barcode"] == precancer), step00.sharing_strategy].itertuples(index=False, name=None))
            intersection = len(primary_set & precancer_set)

            if (clinical_data.loc[patient, step00.sharing_columns[0]] is None) or (clinical_data.loc[patient, step00.sharing_columns[0]] < intersection / len(primary_set)):
                clinical_data.loc[patient, step00.sharing_columns[0]] = intersection / len(primary_set)
                clinical_data.loc[patient, step00.sharing_columns[0] + "-sample"] = precancer

            if (clinical_data.loc[patient, step00.sharing_columns[2]] is None) or (clinical_data.loc[patient, step00.sharing_columns[2]] < intersection / len(precancer_set | primary_set)):
                clinical_data.loc[patient, step00.sharing_columns[2]] = intersection / len(precancer_set | primary_set)
                clinical_data.loc[patient, step00.sharing_columns[2] + "-sample"] = precancer

            if (clinical_data.loc[patient, step00.sharing_columns[4]] is None) or (clinical_data.loc[patient, step00.sharing_columns[4]] < intersection):
                clinical_data.loc[patient, step00.sharing_columns[4]] = intersection
                clinical_data.loc[patient, step00.sharing_columns[4] + "-sample"] = precancer

            if (clinical_data.loc[patient, step00.sharing_columns[6]] is None) or (clinical_data.loc[patient, step00.sharing_columns[6]] < intersection / (len(primary_set) * step00.big / step00.WES_length)):
                clinical_data.loc[patient, step00.sharing_columns[6]] = intersection / (len(primary_set) * step00.big / step00.WES_length)
                clinical_data.loc[patient, step00.sharing_columns[6] + "-sample"] = precancer

    for patient in tqdm.tqdm(patients):
        patient_data = mutect_data.loc[(mutect_data["Patient"] == patient)]
        precancer_list = list(mutect_data.loc[(mutect_data["Patient"] == patient) & ~(mutect_data["Stage"].isin({"Primary"})), "Tumor_Sample_Barcode"])
        primary_set = set(patient_data.loc[(patient_data["Stage"] == "Primary"), step00.sharing_strategy].itertuples(index=False, name=None))

        if (not precancer_list) or (not primary_set):
            continue

        for precancer in precancer_list:
            precancer_set = set(patient_data.loc[(patient_data["Tumor_Sample_Barcode"] == precancer), step00.sharing_strategy].itertuples(index=False, name=None))
            intersection = len(primary_set & precancer_set)

        if (clinical_data.loc[patient, step00.sharing_columns[1]] is None) or (clinical_data.loc[patient, step00.sharing_columns[1]] < intersection / len(primary_set)):
            clinical_data.loc[patient, step00.sharing_columns[1]] = len(primary_set & precancer_set) / len(primary_set)
            clinical_data.loc[patient, step00.sharing_columns[1] + "-sample"] = precancer

        if (clinical_data.loc[patient, step00.sharing_columns[3]] is None) or (clinical_data.loc[patient, step00.sharing_columns[3]] < intersection / len(precancer_set | primary_set)):
            clinical_data.loc[patient, step00.sharing_columns[3]] = intersection / len(precancer_set | primary_set)
            clinical_data.loc[patient, step00.sharing_columns[3] + "-sample"] = precancer

        if (clinical_data.loc[patient, step00.sharing_columns[5]] is None) or (clinical_data.loc[patient, step00.sharing_columns[5]] < intersection):
            clinical_data.loc[patient, step00.sharing_columns[5]] = intersection
            clinical_data.loc[patient, step00.sharing_columns[5] + "-sample"] = precancer

        if (clinical_data.loc[patient, step00.sharing_columns[7]] is None) or (clinical_data.loc[patient, step00.sharing_columns[7]] < intersection / (len(primary_set) * step00.big / step00.WES_length)):
            clinical_data.loc[patient, step00.sharing_columns[7]] = intersection / (len(primary_set) * step00.big / step00.WES_length)
            clinical_data.loc[patient, step00.sharing_columns[7] + "-sample"] = precancer

    clinical_data.dropna(subset=step00.sharing_columns, inplace=True)
    for column in tqdm.tqdm(step00.sharing_columns):
        clinical_data[column] = list(map(float, clinical_data[column]))
    print(clinical_data)

    for column in step00.sharing_columns:
        minimum = numpy.min(clinical_data.loc[(clinical_data["Histology"] == "SQC"), column])
        lower_bound = numpy.quantile(clinical_data.loc[(clinical_data["Histology"] == "SQC"), column], 0.1)
        mean = numpy.mean(clinical_data.loc[(clinical_data["Histology"] == "SQC"), column])
        median = numpy.median(clinical_data.loc[(clinical_data["Histology"] == "SQC"), column])
        higher_bound = numpy.quantile(clinical_data.loc[(clinical_data["Histology"] == "SQC"), column], 0.9)
        maximum = numpy.max(clinical_data.loc[(clinical_data["Histology"] == "SQC"), column])

        print(column, ":", f"min={minimum:.3f}", f"lower={lower_bound:.3f}", f"mean={mean:.3f}", f"median={median:.3f}", f"higher={higher_bound:.3f}", f"max={maximum:.3f}")

    clinical_data.to_csv(args.output, sep="\t")
