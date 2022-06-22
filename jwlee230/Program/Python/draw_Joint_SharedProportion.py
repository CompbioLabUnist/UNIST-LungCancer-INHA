"""
draw_Joint_SharedProportion.py: draw Joint plot with Shared mutation proportion
"""
import argparse
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import scipy.stats
import seaborn
import pandas
import tqdm
import step00

survival_columns = ["Recurrence-Free Survival", "Overall Survival"]


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
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    clinical_data: pandas.DataFrame = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(len(patients))

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))
    args.input.sort(key=step00.sorting)
    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
        mutect_data["Tumor_Sample_Barcode"] = pool.map(step00.get_id, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Patient"] = pool.map(step00.get_patient, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Stage"] = pool.map(step00.get_long_sample_type, mutect_data["Tumor_Sample_Barcode"])
        mutect_data = mutect_data.loc[(mutect_data[step00.nonsynonymous_column].isin(step00.nonsynonymous_mutations))]
    print(mutect_data)

    patients &= set(mutect_data["Patient"])
    print(len(patients))

    clinical_data["Shared Proportion"] = None
    for patient in tqdm.tqdm(patients):
        patient_data = mutect_data.loc[mutect_data["Patient"] == patient]

        stage_set = list(filter(lambda x: x in set(patient_data["Stage"]), step00.long_sample_type_list))

        if ("Primary" not in stage_set) and (len(stage_set) < 2):
            continue

        primary_set = set(patient_data.loc[patient_data["Stage"] == "Primary", step00.sharing_strategy].itertuples(index=False, name=None))
        precancer_set = set(patient_data.loc[patient_data["Stage"] == stage_set[-2], step00.sharing_strategy].itertuples(index=False, name=None))
        clinical_data.loc[patient, "Shared Proportion"] = len(primary_set & precancer_set) / len(primary_set)
    clinical_data.dropna(subset=["Shared Proportion"], inplace=True)
    clinical_data["Shared Proportion"] = list(map(float, clinical_data["Shared Proportion"]))
    print(clinical_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()

    for column in tqdm.tqdm(survival_columns):
        clinical_data[column] = list(map(int, clinical_data[column]))

        r, p = scipy.stats.pearsonr(clinical_data["Shared Proportion"], clinical_data[column])

        g = seaborn.jointplot(data=clinical_data, x="Shared Proportion", y=column, kind="reg", height=24, dropna=True)
        g.fig.text(0.5, 0.75, "r={0:.3f}, p={1:.3f}".format(r, p), color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")

        figures.append(f"{column}.pdf")
        g.savefig(figures[-1])

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
