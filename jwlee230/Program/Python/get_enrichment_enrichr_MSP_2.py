"""
get_enrichment_enrichr_MSP_2.py: get enrichment information with Enrichr between DEG & MSP, but Precancer only
"""
import argparse
import json
import tarfile
import typing
import time
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import requests
import tqdm
import step00


def get_response(url, payload):
    while True:
        if payload is not None:
            response = requests.post(url, files=payload)
        else:
            response = requests.get(url)

        if not response.ok:
            print(f"Response is not ok!! {response.status_code}!!")
            time.sleep(5)
        else:
            data = json.loads(response.text)
            return data


def run(file_name: str, genes: typing.List[str], color: str) -> typing.List[str]:
    tmp_files: typing.List[str] = list()

    if genes:
        gene_set_data = get_response(step00.pathway_addlist_url, {"list": (None, "\n".join(genes)), "description": (None, args.output)})

        raw_data = get_response(f"{step00.pathway_enrichment_url}?userListId={gene_set_data['userListId']}&backgroundType={args.DB}", None)
        enrichment_data = pandas.DataFrame(raw_data[args.DB], columns=step00.pathway_wanted_columns)
        enrichment_data = enrichment_data.loc[(enrichment_data["Adjusted p-value"] < args.padj)]
        enrichment_data["Overlapping genes..."] = list(map(lambda x: ",".join(x) if (len(x) < 10 + 1) else (",".join(x[:10] + ["..."]) + f"({len(x)})"), enrichment_data["Overlapping genes"]))
        enrichment_data["Overlapping genes"] = list(map(lambda x: ",".join(x), enrichment_data["Overlapping genes"]))
        enrichment_data["-log10(P)"] = -1 * numpy.log10(enrichment_data["P-value"])
        enrichment_data["-log10(Padj)"] = -1 * numpy.log10(enrichment_data["Adjusted p-value"])
        enrichment_data["Gene count"] = list(map(lambda x: len(x.split(",")), list(enrichment_data["Overlapping genes"])))

        enrichment_data = enrichment_data.loc[~(enrichment_data["Term name"].str.endswith("disease"))]
        enrichment_data["Rank"] = list(range(len(enrichment_data)))
    else:
        enrichment_data = pandas.DataFrame(columns=step00.pathway_wanted_columns)

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    if (enrichment_data.empty) or (len(enrichment_data) < 1):
        matplotlib.pyplot.text(0.5, 0.5, "Nothing to show...", fontsize="xx-large", color="k", horizontalalignment="center", verticalalignment="center")
        matplotlib.pyplot.xticks([])
        matplotlib.pyplot.yticks([])
    else:
        seaborn.scatterplot(data=enrichment_data, x="-log10(Padj)", y="Rank", size="Gene count", sizes=(100, 1000), hue="Z-score", palette="Reds", legend="brief", edgecolor="black")
        matplotlib.pyplot.axvline(x=-numpy.log10(0.05), color="black", linestyle="--")

        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.yticks(enrichment_data["Rank"], enrichment_data["Term name"], fontsize="xx-small")
        matplotlib.pyplot.xlabel("-log10(Padj)")
        matplotlib.pyplot.ylabel(f"{enrichment_data.shape[0]} pathways")
        matplotlib.pyplot.legend(loc="lower right")
        ax.invert_yaxis()

    matplotlib.pyplot.tight_layout()
    tmp_files.append(f"{file_name}.pdf")
    fig.savefig(tmp_files[-1])
    matplotlib.pyplot.close(fig)

    fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

    if (enrichment_data.empty) or (len(enrichment_data) < 1):
        matplotlib.pyplot.text(0.5, 0.5, "Nothing to show...", fontsize="xx-large", color="k", horizontalalignment="center", verticalalignment="center")
        matplotlib.pyplot.xticks([])
        matplotlib.pyplot.yticks([])
    else:
        matplotlib.pyplot.barh(y=enrichment_data["Rank"], width=enrichment_data["-log10(Padj)"], color=color)
        matplotlib.pyplot.axvline(x=-numpy.log10(0.05), color="black", linestyle="--")

        for index, row in enrichment_data.iterrows():
            matplotlib.pyplot.text(0, y=row["Rank"], s=row["Overlapping genes..."], fontsize="xx-small", horizontalalignment="left", verticalalignment="center")

        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.yticks(enrichment_data["Rank"], enrichment_data["Term name"], fontsize="xx-small")
        matplotlib.pyplot.xlabel("-log10(Padj)")
        matplotlib.pyplot.ylabel(f"{enrichment_data.shape[0]} pathways")
        ax.invert_yaxis()

    matplotlib.pyplot.tight_layout()
    tmp_files.append(f"{file_name}_bar.pdf")
    fig.savefig(tmp_files[-1])
    matplotlib.pyplot.close(fig)

    tmp_files.append(f"{file_name}.tsv")
    enrichment_data.loc[:, step00.pathway_wanted_columns].to_csv(tmp_files[-1], sep="\t", index=False)

    rows = enrichment_data.shape[0]
    enrichment_data = enrichment_data.iloc[:3, :].loc[:, ["Term name", "Overlapping genes...", "Adjusted p-value"]]
    if rows > 3:
        enrichment_data.columns = [f"Term name ({rows})", "Overlapping genes...", "Adjusted p-value"]

    tmp_files.append(f"{file_name}.tex")
    enrichment_data.to_latex(tmp_files[-1], index=False, float_format="%.2e")

    return tmp_files


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input DEG-MSP TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--DB", help="Database name", choices=step00.pathway_gene_set_library, required=True)
    parser.add_argument("--r", help="r-value threshold", type=float, default=0.4)
    parser.add_argument("--slope", help="Slope threshold", type=float, default=7.5)
    parser.add_argument("--padj", help="P-value threshold", type=float, default=0.05)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif not (0 < args.r < 1):
        raise ValueError("r-value must be in (0, 1)!!")
    elif args.slope <= 0:
        raise ValueError("Slope must be positive!!")
    elif not (0 < args.padj < 1):
        raise ValueError("Padj must be (0, 1)!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)

    precancer_stage = "Precancer"
    primary_stage = "Primary"

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    files = list()
    for MSP in tqdm.tqdm(step00.sharing_columns[:2]):
        genes = sorted(set(input_data.loc[(input_data[f"{precancer_stage}-{MSP}-r"] > args.r) & (input_data[f"{precancer_stage}-{MSP}-slope"] > args.slope)].index) - set(input_data.loc[(input_data[f"{primary_stage}-{MSP}-r"] > args.r) & (input_data[f"{primary_stage}-{MSP}-slope"] > args.slope)].index))
        files += run(f"{precancer_stage}-{MSP}-MT-Up", genes, "tab:pink")

        genes = sorted(set(input_data.loc[(input_data[f"{precancer_stage}-{MSP}-r"] < (-1 * args.r)) & (input_data[f"{precancer_stage}-{MSP}-slope"] > args.slope)].index) - set(input_data.loc[(input_data[f"{primary_stage}-{MSP}-r"] < (-1 * args.r)) & (input_data[f"{primary_stage}-{MSP}-slope"] > args.slope)].index))
        files += run(f"{precancer_stage}-{MSP}-MT-Down", genes, "tab:cyan")

    input_data = input_data.loc[list(filter(lambda x: not x.startswith("MT-"), list(input_data.index)))]

    for MSP in tqdm.tqdm(step00.sharing_columns[:2]):
        genes = sorted(set(input_data.loc[(input_data[f"{precancer_stage}-{MSP}-r"] > args.r) & (input_data[f"{precancer_stage}-{MSP}-slope"] > args.slope)].index) - set(input_data.loc[(input_data[f"{primary_stage}-{MSP}-r"] > args.r) & (input_data[f"{primary_stage}-{MSP}-slope"] > args.slope)].index))
        files += run(f"{precancer_stage}-{MSP}-noMT-Up", genes, "tab:pink")

        genes = sorted(set(input_data.loc[(input_data[f"{precancer_stage}-{MSP}-r"] < (-1 * args.r)) & (input_data[f"{precancer_stage}-{MSP}-slope"] > args.slope)].index) - set(input_data.loc[(input_data[f"{primary_stage}-{MSP}-r"] < (-1 * args.r)) & (input_data[f"{primary_stage}-{MSP}-slope"] > args.slope)].index))
        files += run(f"{precancer_stage}-{MSP}-noMT-Down", genes, "tab:cyan")

    with tarfile.open(args.output, "w") as tar:
        for file in tqdm.tqdm(files):
            tar.add(file, arcname=file)
