"""
calculate_beta_PathSeq.py: Calculate beta-diversity for PathSeq
"""
import argparse
import pandas
import skbio.diversity
import skbio.tree

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PathSeq results TSV file", type=str)
    parser.add_argument("tree", help="PathSeq tree NWK file", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--beta", help="Beta-diversity", choices=skbio.diversity.get_beta_diversity_metrics(), required=True)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with TSV!!")
    elif not args.tree.endswith(".nwk"):
        raise ValueError("Tree file must end with NWK!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output file must end with TSV!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0).iloc[:, :-1]
    print(input_data)

    tree = skbio.tree.TreeNode.read(args.tree)

    output_data = skbio.diversity.beta_diversity(args.beta, input_data.to_numpy(), list(input_data.index), tree=tree, otu_ids=list(input_data.columns)).to_data_frame()
    print(output_data)
    output_data.to_csv(args.output, sep="\t")
