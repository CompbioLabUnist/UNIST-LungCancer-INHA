"""
draw_distribution_SharedProportion.py: draw distribution of mutation shared proportion
"""
import argparse
import multiprocessing
import matplotlib
import matplotlib.pyplot
import seaborn
import pandas
import tqdm
import step00

wanted_columns = ["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    clinical_data: pandas.DataFrame = step00.get_clinical_data(args.clinical)
    clinical_data = clinical_data.loc[~(clinical_data["Volume_Doubling_Time"].isna())]
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
    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
        mutect_data["Tumor_Sample_Barcode"] = pool.map(step00.get_id, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Patient"] = pool.map(step00.get_patient, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Stage"] = pool.map(step00.get_long_sample_type, mutect_data["Tumor_Sample_Barcode"])
    print(mutect_data)

    patients &= set(mutect_data["Patient"])

    clinical_data["Shared Proportion"] = 0.0
    for patient in tqdm.tqdm(patients):
        patient_data = mutect_data.loc[mutect_data["Patient"] == patient]

        stage_set = set(patient_data["Stage"])
        assert "Primary" in stage_set
        primary_set = set(patient_data.loc[patient_data["Stage"] == "Primary", wanted_columns].itertuples(index=False, name=None))

        for stage in stage_set:
            if stage == "Primary":
                continue

            precancer_set = set(patient_data.loc[patient_data["Stage"] == stage, wanted_columns].itertuples(index=False, name=None))
            proportion = len(primary_set & precancer_set) / len(primary_set)

            clinical_data.loc[patient, "Shared Proportion"] = max(clinical_data.loc[patient, "Shared Proportion"], proportion)
    print(clinical_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    g = seaborn.displot(clinical_data, x="Shared Proportion", kind="hist", stat="probability", kde=True, height=18, aspect=16 / 9)

    g.tight_layout()
    g.savefig(args.output)
