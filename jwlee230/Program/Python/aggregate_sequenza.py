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
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".txt"), args.input)):
        raise ValueError("INPUT must end with .txt!!")
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

    input_data.columns = list(map(lambda x: x.split("/")[-2], args.input))
    input_data = input_data.T
    input_data["type"] = list(map(lambda x: step00.long_sample_type_dict[x], list(map(step00.get_sample_type, list(input_data.index)))))
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    order = list(filter(lambda x: x in set(input_data["type"]), step00.long_sample_type_list))
    seaborn.scatterplot(data=input_data, x="cellularity", y="ploidy", hue="type", style="type", legend="full", hue_order=order, style_order=order, s=1000, ax=ax)

    matplotlib.pyplot.xlim(-0.1, 1.1)
    matplotlib.pyplot.axvline(x=1, color="k", linestyle="--")
    matplotlib.pyplot.axhline(y=2, color="k", linestyle="--")
    matplotlib.pyplot.xlabel("Cellularity")
    matplotlib.pyplot.ylabel("Ploidy")

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
