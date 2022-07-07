"""
aggregate_sequenza.py: aggregate sequenza results
"""
import argparse
import multiprocessing
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import step00


def get_data(file_name: str) -> pandas.DataFrame:
    data = pandas.read_csv(file_name, header=0, sep="\t")
    data.sort_values("SLPP", ascending=False, inplace=True, ignore_index=True)
    return data.iloc[0, :]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Sequenza output alternative_solutions.txt file(s)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data CSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".txt"), args.input)):
        raise ValueError("INPUT must end with .TXT!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(patients)

    args.input = list(filter(lambda x: step00.get_patient(x.split("/")[-2]) in patients, args.input))

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(get_data, args.input), axis="columns", copy=False)
    print(input_data)

    input_data.columns = list(map(lambda x: x.split("/")[-2], args.input))
    input_data = input_data.T
    input_data["type"] = list(map(lambda x: step00.long_sample_type_dict[x], list(map(step00.get_sample_type, list(input_data.index)))))
    print(input_data)

    sample_list = list(map(lambda x: x.split("/")[-2], args.input))
    stage_set = set(map(step00.get_long_sample_type, sample_list))
    stage_list = list(filter(lambda x: x in stage_set, step00.long_sample_type_list))
    drawing_palette = list(map(lambda x: step00.stage_color_code[x], stage_list))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    order = list(filter(lambda x: x in set(input_data["type"]), step00.long_sample_type_list))
    g = seaborn.JointGrid(data=input_data, x="cellularity", y="ploidy", hue="type", hue_order=order, xlim=(-0.1, 1.1), height=24, ratio=5, palette=drawing_palette)
    g.plot_joint(seaborn.scatterplot, s=1000, legend="full")
    g.plot_marginals(seaborn.histplot, kde=True, stat="probability", multiple="stack")

    g.savefig(args.output)
