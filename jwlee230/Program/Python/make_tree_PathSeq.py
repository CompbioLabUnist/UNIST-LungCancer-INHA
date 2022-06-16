"""
make_tree_PathsSeq.py: make tree from Pathseq data
"""
import argparse
import typing
import pandas
import skbio.tree
import tqdm

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PathSeq results TSV file", type=str, nargs="+")
    parser.add_argument("output", help="Output NWK file", type=str)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("Input file must end with .TSV!!")
    elif not args.output.endswith(".nwk"):
        raise ValueError("Output file must end with .NWK!!")

    lineages: typing.Dict[str, typing.List[str]] = dict()
    for input_file in tqdm.tqdm(args.input):
        data = pandas.read_csv(input_file, sep="\t", dtype=str)
        data["taxonomy"] = list(map(lambda x: x.replace("_", " ").split("|"), data["taxonomy"]))
        data["id"] = list(map(lambda x: x[-1].replace("_", " "), data["taxonomy"]))
        lineages |= dict(data[["id", "taxonomy"]].itertuples(index=False, name=None))

    tree = skbio.tree.TreeNode.from_taxonomy(lineages.items())
    tree.assign_ids()
    for node in tqdm.tqdm(tree.postorder(include_self=True)):
        if node.length is None:
            node.length = 0.0
    tree.write(args.output)
    print(tree.ascii_art())
