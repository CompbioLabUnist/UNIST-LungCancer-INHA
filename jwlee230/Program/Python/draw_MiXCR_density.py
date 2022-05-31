"""
draw_MiXCR_density.py: Draw T-cell density from MiXCR results
"""
import argparse
import itertools
import multiprocessing
import matplotlib
import matplotlib.pyplot
import pandas
import scipy.stats
import seaborn
import statannotations.Annotator
import tqdm
import step00


def read_result(filename: str) -> pandas.DataFrame:
    data = pandas.read_csv(filename, sep="\t")
    data = data.loc[(data["allVHitsWithScore"].str.contains("TRB"))]
    data["sample"] = step00.get_id(filename)
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="MiXCR result TXT file", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinical data data CSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".txt"), args.input)):
        raise ValueError("Input must end with .TXT!!")
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

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(pool.map(read_result, args.input), join="outer", ignore_index=True)
    print(input_data)

    output_data = pandas.DataFrame(index=list(map(step00.get_id, args.input)))
    output_data["Stage"] = list(map(step00.get_long_sample_type, list(output_data.index)))
    output_data["Density"] = 0
    for index, row in tqdm.tqdm(input_data.iterrows()):
        output_data.loc[row["sample"], "Density"] += row["cloneFraction"]
    print(output_data)

    order = list(filter(lambda x: x in list(output_data["Stage"]), step00.long_sample_type_list))
    print(order)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    try:
        stat, p = scipy.stats.kruskal(*[output_data.loc[(output_data["Stage"] == stage), "Density"] for stage in order])
    except ValueError:
        _, p = 0.0, 1.0

    seaborn.violinplot(data=output_data, x="Stage", y="Density", order=order, palette=step00.stage_color_code, cut=1, linewidth=5, ax=ax)
    statannotations.Annotator.Annotator(ax, list(itertools.combinations(order, r=2)), data=output_data, x="Stage", y="Density", order=order, palette=step00.stage_color_code).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.title(f"Kruskal-Wallis p={p:.3f}")
    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
