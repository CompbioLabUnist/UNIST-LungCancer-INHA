"""
draw_KM_SharedProportion.py: Draw Kaplan-Meier plot with mutation shared proportion
"""
import argparse
import multiprocessing
import tarfile
import lifelines
import lifelines.statistics
import matplotlib
import matplotlib.pyplot
import numpy
import seaborn
import pandas
import tqdm
import step00

wanted_columns = ["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]
survival_columns = ["Recurrence-Free Survival", "Overall Survival"]
threshold = 365 * 5


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)
    parser.add_argument("--cutting", help="Cutting follow-up up to 5 year", action="store_true", default=False)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    group_divide = parser.add_mutually_exclusive_group(required=True)
    group_divide.add_argument("--median", help="Divide as median", action="store_true", default=False)
    group_divide.add_argument("--mean", help="Divide as mean", action="store_true", default=False)

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
    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
        mutect_data["Tumor_Sample_Barcode"] = pool.map(step00.get_id, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Patient"] = pool.map(step00.get_patient, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Stage"] = pool.map(step00.get_long_sample_type, mutect_data["Tumor_Sample_Barcode"])
    print(mutect_data)

    patients &= set(mutect_data["Patient"])
    print(len(patients))

    clinical_data["Shared Proportion"] = 0.0
    for patient in tqdm.tqdm(patients):
        patient_data = mutect_data.loc[mutect_data["Patient"] == patient]

        stage_set = set(patient_data["Stage"])
        if "Primary" not in stage_set:
            continue
        primary_set = set(patient_data.loc[patient_data["Stage"] == "Primary", wanted_columns].itertuples(index=False, name=None))

        for stage in stage_set:
            if stage == "Primary":
                continue

            precancer_set = set(patient_data.loc[patient_data["Stage"] == stage, wanted_columns].itertuples(index=False, name=None))
            proportion = len(primary_set & precancer_set) / len(primary_set)

            clinical_data.loc[patient, "Shared Proportion"] = max(clinical_data.loc[patient, "Shared Proportion"], proportion)

    if args.cutting:
        for column in tqdm.tqdm(survival_columns):
            clinical_data[column] = list(map(lambda x: threshold if (x > threshold) else x, clinical_data[column]))
    print(clinical_data)

    if args.median:
        median = numpy.median(clinical_data["Shared Proportion"])
        lower_data = clinical_data.loc[clinical_data["Shared Proportion"] < median]
        higher_data = clinical_data.loc[clinical_data["Shared Proportion"] >= median]
    elif args.mean:
        mean = numpy.mean(clinical_data["Shared Proportion"])
        lower_data = clinical_data.loc[clinical_data["Shared Proportion"] < mean]
        higher_data = clinical_data.loc[clinical_data["Shared Proportion"] >= mean]
    else:
        raise Exception("Something went wrong!!")
    print(lower_data)
    print(higher_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()

    for column in tqdm.tqdm(survival_columns):
        fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

        kmf = lifelines.KaplanMeierFitter()

        kmf.fit(lower_data[column], label=f"Lower Shared Proportion ({len(lower_data)} patients)")
        kmf.plot(ax=ax, ci_show=False, c="tab:blue")

        kmf.fit(higher_data[column], label=f"Higher Shared Proportion ({len(higher_data)} patients)")
        kmf.plot(ax=ax, ci_show=False, c="tab:red")

        p_value = lifelines.statistics.logrank_test(lower_data[column], higher_data[column]).p_value
        matplotlib.pyplot.text(max(clinical_data[column]) / 2, 0.75, f"p={p_value:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")

        matplotlib.pyplot.xlabel(f"{column} (Days)")
        matplotlib.pyplot.ylabel("Survival Rate")
        matplotlib.pyplot.tight_layout()

        figures.append(f"{column.replace(' ', '-')}.pdf")
        matplotlib.pyplot.savefig(figures[-1])
        matplotlib.pyplot.close()

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
