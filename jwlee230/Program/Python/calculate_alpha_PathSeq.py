"""
calculate_alpha_PathSeq.py: Calculate alpha-diversity for PathSeq
"""
import argparse
import pandas
import skbio.diversity
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="PathSeq results TSV file", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with TSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output file must end with TSV!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0).iloc[:, :-1]
    print(input_data)

    output_data = pandas.DataFrame(index=input_data.index)

    for alpha in tqdm.tqdm(skbio.diversity.get_alpha_diversity_metrics()):
        if alpha.endswith("_ci"):
            continue
        try:
            output_data[alpha] = skbio.diversity.alpha_diversity(alpha, input_data.to_numpy(), list(input_data.index))
        except Exception:
            continue

    output_data["Subtype"] = list(map(step00.get_long_sample_type, list(output_data.index)))
    print(output_data)
    output_data.to_csv(args.output, sep="\t")
