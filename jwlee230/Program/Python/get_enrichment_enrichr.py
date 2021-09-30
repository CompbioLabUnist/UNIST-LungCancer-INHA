"""
get_enrichment_enrichr.py: get enrichment information with Enrichr
"""
import argparse
import json
import numpy
import pandas
import requests

addlist_url = "https://maayanlab.cloud/Enrichr/addList"
enrichment_url = "https://maayanlab.cloud/Enrichr/enrich"
gene_set_library = ["KEGG_2021_Human", "MSigDB_Oncogenic_Signatures"]


def get_response(url, payload):
    if payload is not None:
        response = requests.post(url, files=payload)
    else:
        response = requests.get(url)

    if not response.ok:
        raise Exception("Response is not Ok!!")

    data = json.loads(response.text)
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("DEG", help="DEG TSV file", type=str)
    parser.add_argument("output", help="Output basename file", type=str)
    parser.add_argument("--DB", help="Database name", choices=gene_set_library, required=True)
    parser.add_argument("--padj", help="P-value threshold", type=float, default=0.05)
    parser.add_argument("--fold", help="Fold change threshold", type=float, default=2)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--up", help="Test up-regulated pathway", action="store_true", default=False)
    group.add_argument("--down", help="Test down-regulated pathway", action="store_true", default=False)

    args = parser.parse_args()

    if not args.DEG.endswith(".tsv"):
        raise ValueError("DEG must end with .TSV!!")

    DEG_data = pandas.read_csv(args.DEG, sep="\t", header=0, names=["gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"], index_col="gene_id").dropna(axis="index", how="any")
    DEG_data["-log(Padj)"] = -1 * numpy.log10(DEG_data["padj"], dtype=float)
    print(DEG_data)

    if args.up:
        genes = list(DEG_data.loc[(DEG_data["log2FoldChange"] >= numpy.log2(args.fold)) & (DEG_data["padj"] < args.padj) & (DEG_data["pvalue"] < args.padj), ["log2FoldChange", "-log(Padj)"]].index)
    elif args.down:
        genes = list(DEG_data.loc[(DEG_data["log2FoldChange"] <= -1 * numpy.log2(args.fold)) & (DEG_data["padj"] < args.padj) & (DEG_data["pvalue"] < args.padj), ["log2FoldChange", "-log(Padj)"]].index)
    else:
        raise Exception("Something went wrong!!")
    print(genes)

    if genes:
        gene_set_data = get_response(addlist_url, {"list": (None, "\n".join(genes)), "description": (None, args.output)})
        print(gene_set_data)

        raw_data = get_response("{0}?userListId={1}&backgroundType={2}".format(enrichment_url, gene_set_data["userListId"], args.DB), None)
        enrichment_data = pandas.DataFrame(raw_data[args.DB], columns=["Rank", "Term name", "P-value", "Z-score", "Combined score", "Overlapping genes", "Adjusted p-value", "Old p-value", "Old adjusted p-value"])
        enrichment_data["Overlapping genes"] = list(map(lambda x: ",".join(sorted(x)), enrichment_data["Overlapping genes"]))
        enrichment_data = enrichment_data.loc[(enrichment_data["P-value"] < args.padj) & (enrichment_data["Adjusted p-value"] < args.padj)]
        print(enrichment_data)
    else:
        enrichment_data = pandas.DataFrame(columns=["Rank", "Term name", "P-value", "Z-score", "Combined score", "Overlapping genes", "Adjusted p-value", "Old p-value", "Old adjusted p-value"])

    enrichment_data.to_csv(args.output + ".tsv", sep="\t", index=False)
    if enrichment_data.iloc[:3, :].empty:
        pandas.DataFrame(columns=["Term name", "Adjusted p-value"]).to_latex(args.output + ".tex", index=False, float_format="%.2e")
    else:
        enrichment_data.iloc[:3, :].loc[:, ["Term name", "Adjusted p-value"]].to_latex(args.output + ".tex", index=False, float_format="%.2e")
