"""
compare_MuationSharedProportion_Sequenza.py: compare Mutation Shared Proportion vs. Sequenza
"""
import argparse
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import scipy.stats
import seaborn
import tqdm
import step00


def get_data(file_name: str) -> pandas.DataFrame:
    data = pandas.read_csv(file_name, header=0, sep="\t")
    data.sort_values("SLPP", ascending=False, inplace=True, ignore_index=True)
    return data.iloc[0, :]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Sequenza output alternative_solutions.txt file(s)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data with Mutaion Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".txt"), args.input)):
        raise ValueError("INPUT must end with .TXT!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
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
    input_data["Stage"] = list(map(lambda x: step00.long_sample_type_dict[x], list(map(step00.get_sample_type, list(input_data.index)))))
    input_data["Patient"] = list(map(step00.get_patient, list(input_data.index)))
    input_data["Mutation Shared Proportion"] = list(map(lambda x: clinical_data.loc[x, "Shared Proportion"], list(input_data["Patient"])))
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    sample_list = list(map(lambda x: x.split("/")[-2], args.input))
    stage_set = set(map(step00.get_long_sample_type, sample_list))
    stage_list = list(filter(lambda x: x in stage_set, step00.long_sample_type_list))

    figures = list()
    for stage in tqdm.tqdm(stage_list):
        drawing_data = input_data.loc[(input_data["Stage"] == stage)]

        r, p = scipy.stats.pearsonr(drawing_data["Mutation Shared Proportion"], drawing_data["cellularity"])

        g = seaborn.jointplot(data=drawing_data, x="Mutation Shared Proportion", y="cellularity", kind="reg", height=24, ratio=5, color=step00.stage_color_code[stage])
        g.fig.text(0.5, 0.75, "r={0:.3f}, p={1:.3f}".format(r, p), color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"}, fontfamily="monospace")

        figures.append(f"{stage}.pdf")
        g.savefig(figures[-1])

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
