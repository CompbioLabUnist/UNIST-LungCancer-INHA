"""
draw_KM_SharedProportion.py: Draw Kaplan-Meier plot with mutation shared proportion
"""
import argparse
import itertools
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

survival_columns = ["Recurrence-Free Survival", "Overall Survival"]
threshold = 365 * 5


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutation Sharing Proportion input TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cutting", help="Cutting follow-up up to 5 year", action="store_true", default=False)
    parser.add_argument("--percentage", help="Percentage threshold", type=float, default=0.33)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .tsv!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    clinical_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "SQC")]
    elif args.ADC:
        clinical_data = clinical_data.loc[(clinical_data["Histology"] == "ADC")]
    else:
        raise Exception("Something went wrong!!")
    patients = set(clinical_data.index)
    print(len(patients))

    if args.cutting:
        for column in tqdm.tqdm(survival_columns):
            clinical_data[column] = list(map(lambda x: threshold if (x > threshold) else x, clinical_data[column]))
    print(clinical_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()

    for MSP, column in tqdm.tqdm(list(itertools.product(step00.sharing_columns[:2], survival_columns))):
        precancer_list = sorted(clinical_data[f"{MSP}-sample"], key=lambda x: clinical_data.loc[step00.get_patient(x), MSP])
        primary_list = list(map(step00.get_paired_primary, precancer_list))
        patient_list = list(map(step00.get_patient, precancer_list))

        MSP_Q1 = numpy.quantile(clinical_data[MSP], args.percentage)
        MSP_Q3 = numpy.quantile(clinical_data[MSP], 1.0 - args.percentage)

        MSP_Q_list = list(map(lambda x: "PSM-L" if (clinical_data.loc[x, MSP] <= MSP_Q1) else ("PSM-H" if (clinical_data.loc[x, MSP] >= MSP_Q3) else "None"), patient_list))
        MSP_L_list = list(map(lambda x: x[1], list(filter(lambda x: x[0] == "PSM-L", zip(MSP_Q_list, precancer_list)))))
        MSP_H_list = list(map(lambda x: x[1], list(filter(lambda x: x[0] == "PSM-H", zip(MSP_Q_list, precancer_list)))))

        if column == "Recurrence-Free Survival":
            lower_events = list(map(lambda x: (x[0] == "1") or (x[1] == "1"), clinical_data.loc[list(map(step00.get_patient, MSP_L_list)), ["Death", "Recurrence"]].itertuples(index=False, name=None)))
            higher_events = list(map(lambda x: (x[0] == "1") or (x[1] == "1"), clinical_data.loc[list(map(step00.get_patient, MSP_H_list)), ["Death", "Recurrence"]].itertuples(index=False, name=None)))
        else:
            lower_events = list(map(lambda x: x == "1", clinical_data.loc[list(map(step00.get_patient, MSP_L_list)), "Death"]))
            higher_events = list(map(lambda x: x == "1", clinical_data.loc[list(map(step00.get_patient, MSP_H_list)), "Death"]))

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        kmf = lifelines.KaplanMeierFitter()

        kmf.fit(clinical_data.loc[list(map(step00.get_patient, MSP_L_list)), column], event_observed=lower_events, label=f"PSM-L ({len(MSP_L_list)} patients)")
        kmf.plot(ax=ax, ci_show=False, c="tab:blue")

        kmf.fit(clinical_data.loc[list(map(step00.get_patient, MSP_H_list)), column], event_observed=higher_events, label=f"PSM-H ({len(MSP_H_list)} patients)")
        kmf.plot(ax=ax, ci_show=False, c="tab:red")

        p_value = lifelines.statistics.logrank_test(clinical_data.loc[list(map(step00.get_patient, MSP_L_list)), column], clinical_data.loc[list(map(step00.get_patient, MSP_H_list)), column], event_observed_A=lower_events, event_observed_B=higher_events).p_value
        matplotlib.pyplot.text(max(clinical_data[column]) / 2, 0.9, f"p={p_value:.3f}", color="black", fontsize="small", horizontalalignment="center", verticalalignment="center")

        matplotlib.pyplot.xlabel(f"{column} (Days)")
        matplotlib.pyplot.ylabel("Survival Rate")
        matplotlib.pyplot.tight_layout()

        figures.append(f"{MSP}-{column}.pdf")
        matplotlib.pyplot.savefig(figures[-1])
        matplotlib.pyplot.close()

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
