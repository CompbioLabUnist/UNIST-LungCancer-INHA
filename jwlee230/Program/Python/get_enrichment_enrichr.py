"""
get_enrichment_enrichr.py: get enrichment information with Enrichr
"""
import argparse
import json
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import requests
import tqdm
import step00

wanted_columns = ["Rank", "Term name", "P-value", "Z-score", "Combined score", "Overlapping genes", "Adjusted p-value", "Old p-value", "Old adjusted p-value"]
addlist_url = "https://maayanlab.cloud/Enrichr/addList"
enrichment_url = "https://maayanlab.cloud/Enrichr/enrich"
gene_set_library = ["KEGG_2021_Human", "MSigDB_Oncogenic_Signatures"]


def get_response(url, payload):
    if payload is not None:
        response = requests.post(url, files=payload)
    else:
        response = requests.get(url)

    if not response.ok:
        raise Exception("Response is not Ok!! {0}!!".format(response.status_code))

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
    elif not (0 < args.padj < 1):
        raise ValueError("Padj must be (0, 1)!!")

    DEG_data = pandas.read_csv(args.DEG, sep="\t", index_col=0).dropna(axis="index", how="any")
    DEG_data["-log(Padj)"] = -1 * numpy.log10(DEG_data["padj"], dtype=float)
    print(DEG_data)

    if args.up:
        genes = list(DEG_data.loc[(DEG_data["log2FoldChange"] >= numpy.log2(args.fold)) & (DEG_data["padj"] < args.padj) & (DEG_data["pvalue"] < args.padj), ["log2FoldChange", "-log(Padj)"]].index)
    elif args.down:
        genes = list(DEG_data.loc[(DEG_data["log2FoldChange"] <= -1 * numpy.log2(args.fold)) & (DEG_data["padj"] < args.padj) & (DEG_data["pvalue"] < args.padj), ["log2FoldChange", "-log(Padj)"]].index)
    else:
        raise Exception("Something went wrong!!")

    print(sorted(genes))

    if genes:
        gene_set_data = get_response(addlist_url, {"list": (None, "\n".join(genes)), "description": (None, args.output)})
        print(gene_set_data)

        raw_data = get_response("{0}?userListId={1}&backgroundType={2}".format(enrichment_url, gene_set_data["userListId"], args.DB), None)
        enrichment_data = pandas.DataFrame(raw_data[args.DB], columns=wanted_columns)
        enrichment_data = enrichment_data.loc[(enrichment_data["P-value"] < args.padj) & (enrichment_data["Adjusted p-value"] < args.padj)]
        enrichment_data["Overlapping genes..."] = list(map(lambda x: ",".join(x) if (len(x) < 4) else (",".join(x[:3] + ["..."]) + "({0})".format(len(x))), enrichment_data["Overlapping genes"]))
        enrichment_data["Overlapping genes"] = list(map(lambda x: ",".join(x), enrichment_data["Overlapping genes"]))
    else:
        enrichment_data = pandas.DataFrame(columns=wanted_columns)
    print(enrichment_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(64, 18))

    if enrichment_data.empty:
        enrichment_data = pandas.DataFrame(columns=wanted_columns + ["Overlapping genes..."], index=[0], data=[["None"] + [""] * (len(wanted_columns))])

        matplotlib.pyplot.text(0.5, 0.5, "Nothing to show...", fontsize=step00.matplotlib_parameters["axes.titlesize"], color="k", horizontalalignment="center", verticalalignment="center")
        matplotlib.pyplot.xticks([])
        matplotlib.pyplot.yticks([])
    else:
        enrichment_data["-log10(P)"] = -1 * numpy.log10(enrichment_data["P-value"])
        enrichment_data["-log10(Padj)"] = -1 * numpy.log10(enrichment_data["Adjusted p-value"])
        print(enrichment_data)

        rows = enrichment_data.shape[0]
        drawing_data = enrichment_data.iloc[:10, :]

        if args.up:
            ax.barh(range(drawing_data.shape[0]), drawing_data["-log10(Padj)"], color="tab:pink")
        elif args.down:
            ax.barh(range(drawing_data.shape[0]), drawing_data["-log10(Padj)"], color="tab:cyan")

        for index, row in tqdm.tqdm(drawing_data.iterrows()):
            matplotlib.pyplot.text(0, index, "{0}: {1}".format(row["Term name"], row["Overlapping genes..."]), color="k", horizontalalignment="left", verticalalignment="center")

        matplotlib.pyplot.yticks([])
        matplotlib.pyplot.xlabel("-log10(Padj)")
        matplotlib.pyplot.ylabel("{0} pathways".format(rows))
        matplotlib.pyplot.ylim(-1, 10)
        matplotlib.pyplot.axvline(x=-1 * numpy.log10(args.padj), linestyle="--", color="black", alpha=0.5)
        ax.invert_yaxis()
        matplotlib.pyplot.tight_layout()

    fig.savefig(args.output + ".pdf")
    matplotlib.pyplot.close(fig)

    enrichment_data.loc[:, wanted_columns].to_csv(args.output + ".tsv", sep="\t", index=False)
    rows = enrichment_data.shape[0]
    enrichment_data = enrichment_data.iloc[:3, :].loc[:, ["Term name", "Overlapping genes...", "Adjusted p-value"]]
    if rows > 3:
        enrichment_data.columns = ["Term name ({0})".format(rows), "Overlapping genes...", "Adjusted p-value"]
    enrichment_data.to_latex(args.output + ".tex", index=False, float_format="%.2e")
