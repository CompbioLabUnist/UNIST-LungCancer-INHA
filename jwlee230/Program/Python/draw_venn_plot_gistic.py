"""
draw_venn_plot_gistic.py: draw venn diagram plot with gistic peak
"""
import argparse
import re
import matplotlib
import matplotlib.pyplot
import pandas
import venn
import tqdm
import step00


def generate_logics(n_sets):
    for i in range(1, 2 ** n_sets):
        yield bin(i)[2:].zfill(n_sets)


def sorting_cytoband(s):
    s_split = list(map(int, re.split(r"p|q|\.", s)))

    if "p" in s:
        s_split.insert(1, "p")
    elif "q" in s:
        s_split.insert(1, "q")
    else:
        raise Exception("Something went wrong!!")

    return s_split


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input Gistic all_lesions.conf_99.txt file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output file basename", type=str)
    parser.add_argument("--annotation", help="Annotation for venn diagram", type=str, nargs="+", required=True)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--amplification", help="Draw Amplification pathway", action="store_true", default=False)
    group.add_argument("--deletion", help="Draw Deletion pathway", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith("/all_lesions.conf_99.txt"), args.input)):
        raise ValueError("Input is not valid!!")
    elif len(args.input) != len(args.annotation):
        raise ValueError("Annotation must be one-to-one upon DEG!!")

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
        input_data[annotation] = set(filter(lambda x: (not x.startswith("X")) and (not x.startswith("Y")), input_data[annotation]))
    print(input_data)

    every_cytoband = sorted(set.union(*list(input_data.values())), key=sorting_cytoband)
    print(every_cytoband)

    if every_cytoband:
        output_data = pandas.DataFrame(data=[["" for x in args.annotation] for y in every_cytoband], index=every_cytoband, columns=args.annotation, dtype=str)
        for annotation in tqdm.tqdm(args.annotation):
            output_data.loc[input_data[annotation], annotation] = "*"
    else:
        output_data = pandas.DataFrame(data=[["" for x in args.annotation]], index=[""], columns=args.annotation, dtype=str)

    print(output_data)

    output_data.to_latex(args.output + ".tex", column_format="l" + "c" * len(args.annotation))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    try:
        venn.venn(input_data, ax=ax, fmt=step00.venn_format, fontsize=step00.matplotlib_parameters["legend.fontsize"], legend_loc="upper left")
    except ZeroDivisionError:
        pass
    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output + ".pdf")
    matplotlib.pyplot.close(fig)

    datasets = list(input_data.values())
    dataset_union = set(every_cytoband)
    for logic in generate_logics(len(args.annotation)):
        included_sets = [datasets[i] for i in range(len(args.annotation)) if logic[i] == "1"]
        excluded_sets = [datasets[i] for i in range(len(args.annotation)) if logic[i] == "0"]

        print("+".join(list(map(lambda y: y[0], list(filter(lambda x: x[1] == "1", zip(args.annotation, list(logic))))))))
        print(sorted((dataset_union & set.intersection(*included_sets)) - set.union(set(), *excluded_sets), key=sorting_cytoband))
        print()
