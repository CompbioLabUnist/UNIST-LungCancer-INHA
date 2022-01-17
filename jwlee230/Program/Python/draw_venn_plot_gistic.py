"""
draw_venn_plot_gistic.py: draw venn diagram plot with gistic peak
"""
import argparse
import matplotlib
import matplotlib.pyplot
import pandas
import venn
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input Gistic all_lesions.conf_99.txt file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--annotation", help="Annotation for venn diagram", type=str, nargs="+", required=True)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--amplification", help="Draw Amplification pathway", action="store_true", default=False)
    group.add_argument("--deletion", help="Draw Deletion pathway", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith("/all_lesions.conf_99.txt"), args.input)):
        raise ValueError("Input is not valid!!")
    if len(args.input) != len(args.annotation):
        raise ValueError("Annotation must be one-to-one upon DEG!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    input_data = dict()
    for annotation, input_file in tqdm.tqdm(list(zip(args.annotation, args.input))):
        data = pandas.read_csv(input_file, sep="\t")
        data["Descriptor"] = list(map(lambda x: x.strip(), data["Descriptor"]))
        if args.amplification:
            input_data[annotation] = set(data.loc[(data["Unique Name"].str.contains("Amplification Peak")) & ~(data["Unique Name"].str.contains("CN")), "Descriptor"])
        elif args.deletion:
            input_data[annotation] = set(data.loc[(data["Unique Name"].str.contains("Deletion Peak")) & ~(data["Unique Name"].str.contains("CN")), "Descriptor"])
        else:
            raise Exception("Something went wrong!!")
        print(data)
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    venn.venn(input_data, ax=ax, fmt=step00.venn_format, fontsize=step00.matplotlib_parameters["legend.fontsize"], legend_loc="upper left")
    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
