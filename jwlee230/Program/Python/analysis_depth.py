"""
analysis_depth.py: analysis samtools depth
"""
import argparse
import multiprocessing
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import numpy
import pandas
import step00


def get_depth(file_name: str) -> float:
    """
    get_depth: get depth of coverage from file
    """
    print(file_name)
    data = pandas.read_csv(file_name, sep="\t", names=["CHROM", "POS", "Depth"], usecols=["Depth"], skiprows=1)
    return numpy.mean(data["Depth"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Samtools depth TSV files", type=str, nargs="+")
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .pdf!!")
    elif args.cpus < 0:
        raise ValueError("CPUS must be positive!!")

    args.input.sort()
    IDs = list(map(lambda x: x.split(".")[0], list(map(lambda x: x.split("/")[-1], args.input))))

    with multiprocessing.Pool(processes=args.cpus) as pool:
        depths = pool.map(get_depth, args.input)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

    colors = list(matplotlib.colors.TABLEAU_COLORS.keys())
    colors_length = len(colors)
    used_color = 1
    previous_patient = ""
    for i, (ID, depth) in enumerate(zip(IDs, depths)):
        patient = step00.get_patient(ID)
        if previous_patient != patient:
            previous_patient = patient
            used_color += 1

        matplotlib.pyplot.bar(i, depth, color=matplotlib.colors.TABLEAU_COLORS[colors[used_color % colors_length]])

    matplotlib.pyplot.xticks([])
    matplotlib.pyplot.xlabel("Total {0} Patients".format(used_color))
    matplotlib.pyplot.ylabel("Depth")
    matplotlib.pyplot.title("Depth of coverage")
    matplotlib.pyplot.grid(True)

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
