"""
aggregate_sequenza.py: aggregate sequenza results
"""
import argparse
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import step00


def untar(file_name: str) -> pandas.DataFrame:
    with tarfile.open(file_name, "r:gz") as tar:
        txt_file = list(filter(lambda x: x.endswith("_alternative_solutions.txt"), tar.getnames()))[0]
        print(txt_file)
        tar.extract(txt_file, path=step00.tmpfs)

    data = pandas.read_csv(step00.tmpfs + "/" + txt_file, header=0, sep="\t")
    data.sort_values("SLPP", ascending=False, inplace=True, ignore_index=True)
    return data.iloc[0, :]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Sequenza tar.gz file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tar.gz"), args.input)):
        raise ValueError("INPUT must end with .tar.gz!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    with multiprocessing.Pool(args.cpus) as pool:
        input_data = pandas.concat(objs=pool.map(untar, args.input), axis="columns", copy=False)

    input_data.columns = list(map(lambda x: x.split("/")[-1].split(".")[0], args.input))
    input_data = input_data.T
    input_data["type"] = list(map(lambda x: step00.long_sample_type_dict[x], list(map(step00.get_sample_type, list(input_data.index)))))

    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
    seaborn.scatterplot(data=input_data, x="cellularity", y="ploidy", hue="type", style="type", legend="full", hue_order=sorted(set(input_data["type"])), style_order=sorted(set(input_data["type"])), s=1000, ax=ax)
    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
