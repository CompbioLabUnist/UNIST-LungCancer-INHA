"""
calculate_MutationSharedProportion.py: calculate Mutation Shared Proportion
"""
import argparse
import multiprocessing
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

    clinical_data[step00.sharing_columns[0]] = None
    for patient in tqdm.tqdm(patients):
        patient_data = mutect_data.loc[(mutect_data["Patient"] == patient) & (mutect_data[step00.nonsynonymous_column].isin(step00.nonsynonymous_mutations))]
        stage_set = list(filter(lambda x: x in set(patient_data["Stage"]), step00.long_sample_type_list))

        if ("Primary" not in stage_set) and (len(stage_set) < 2):
            continue

        primary_set = set(patient_data.loc[patient_data["Stage"] == "Primary", step00.sharing_strategy].itertuples(index=False, name=None))
        precancer_set = set(patient_data.loc[patient_data["Stage"] == stage_set[-2], step00.sharing_strategy].itertuples(index=False, name=None))
        clinical_data.loc[patient, step00.sharing_columns[0]] = len(primary_set & precancer_set) / len(primary_set)

    clinical_data[step00.sharing_columns[1]] = None
    for patient in tqdm.tqdm(patients):
        patient_data = mutect_data.loc[(mutect_data["Patient"] == patient)]
        stage_set = list(filter(lambda x: x in set(patient_data["Stage"]), step00.long_sample_type_list))

        if ("Primary" not in stage_set) and (len(stage_set) < 2):
            continue

        primary_set = set(patient_data.loc[patient_data["Stage"] == "Primary", step00.sharing_strategy].itertuples(index=False, name=None))
        precancer_set = set(patient_data.loc[patient_data["Stage"] == stage_set[-2], step00.sharing_strategy].itertuples(index=False, name=None))
        clinical_data.loc[patient, step00.sharing_columns[1]] = len(primary_set & precancer_set) / len(primary_set)

    clinical_data[step00.sharing_columns[2]] = None
    for patient in tqdm.tqdm(patients):
        patient_data = mutect_data.loc[(mutect_data["Patient"] == patient) & (mutect_data[step00.nonsynonymous_column].isin(step00.nonsynonymous_mutations))]
        stage_set = list(filter(lambda x: x in set(patient_data["Stage"]), step00.long_sample_type_list))

        if ("Primary" not in stage_set) and (len(stage_set) < 2):
            continue

        primary_set = set(patient_data.loc[patient_data["Stage"] == "Primary", step00.sharing_strategy].itertuples(index=False, name=None))
        precancer_set = set(patient_data.loc[patient_data["Stage"] == stage_set[-2], step00.sharing_strategy].itertuples(index=False, name=None))
        clinical_data.loc[patient, step00.sharing_columns[2]] = len(primary_set & precancer_set) / len(primary_set | precancer_set)

    clinical_data[step00.sharing_columns[3]] = None
    for patient in tqdm.tqdm(patients):
        patient_data = mutect_data.loc[(mutect_data["Patient"] == patient)]
        stage_set = list(filter(lambda x: x in set(patient_data["Stage"]), step00.long_sample_type_list))

        if ("Primary" not in stage_set) and (len(stage_set) < 2):
            continue

        primary_set = set(patient_data.loc[patient_data["Stage"] == "Primary", step00.sharing_strategy].itertuples(index=False, name=None))
        precancer_set = set(patient_data.loc[patient_data["Stage"] == stage_set[-2], step00.sharing_strategy].itertuples(index=False, name=None))
        clinical_data.loc[patient, step00.sharing_columns[3]] = len(primary_set & precancer_set) / len(primary_set | precancer_set)

    clinical_data.dropna(subset=step00.sharing_columns, inplace=True)
    for column in tqdm.tqdm(step00.sharing_columns):
        clinical_data[column] = list(map(float, clinical_data[column]))
    print(clinical_data)

    clinical_data.to_csv(args.output, sep="\t")
