"""
convert_mutect2_pyclone.py: convert Mutect2 result to PyClone input format
"""
import argparse
import typing
from adjustText import adjust_text
import gtfparse
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("loci", help="Pyclone loci TSV file", type=str)
    parser.add_argument("census", help="Cancer gene census CSV file", type=str)
    parser.add_argument("known", help="Known gene GTF.gz file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.05)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--CCF", help="Draw CCF (cancer cell fraction) plot", action="store_true", default=False)
    group.add_argument("--VAF", help="Draw VAF (variant allele frequency) plot", action="store_true", default=False)

    args = parser.parse_args()

    if not args.loci.endswith(".tsv"):
        raise ValueError("Loci must end with .TSV!!")
    elif not args.census.endswith(".csv"):
        raise ValueError("Census must end with .CSV!!")
    elif not args.known.endswith(".gtf.gz"):
        raise ValueError("Known must end with .GTF.gz!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif not (0 < args.p < 1):
        raise ValueError("P-value must be (0, 1)!!")

    first_name = step00.get_id(args.output)
    second_name = step00.get_paired_primary(first_name)

    try:
        pyclone_data = pandas.read_csv(args.loci, sep="\t")
    except pandas.errors.EmptyDataError:
        pyclone_data = pandas.DataFrame(columns=["mutation_id", "sample_id", "cluster_id", "cellular_prevalence", "cellular_prevalence_std", "variant_allele_frequency"])

    pyclone_data.sort_values(by="sample_id", inplace=True, ignore_index=True)
    pyclone_data["chromosome"] = list(map(lambda x: x.split(":")[0], pyclone_data["mutation_id"]))
    pyclone_data["start"] = list(map(lambda x: int(x.split(":")[1]), pyclone_data["mutation_id"]))
    pyclone_data["end"] = list(map(lambda x: int(x.split(":")[2]), pyclone_data["mutation_id"]))
    print(pyclone_data)

    census_data = pandas.read_csv(args.census)
    census_set = set(census_data["Gene Symbol"])
    print(census_data)

    gtf_data = gtfparse.read_gtf(args.known)
    print(gtf_data)

    gene_names: typing.List[typing.List[str]] = list()
    for index, row in tqdm.tqdm(pyclone_data.iterrows()):
        tmp_data = gtf_data.loc[(gtf_data["seqname"] == row["chromosome"]) & (gtf_data["start"] <= row["start"]) & (row["end"] <= gtf_data["end"]), "gene_id"]

        if tmp_data.empty:
            gene_names.append([])
        else:
            gene_names.append(sorted(set(tmp_data.to_numpy())))

    pyclone_data["gene"] = gene_names
    pyclone_data = pyclone_data.explode("gene", ignore_index=True)
    print(pyclone_data)

    pyclone_data["gene_census"] = list(map(lambda x: x in census_set, pyclone_data["gene"]))
    pyclone_data["mutation"] = list(map(lambda x: not x.endswith("nan"), pyclone_data["mutation_id"]))
    pyclone_data.rename(columns={"cluster_id": "Cluster ID"}, inplace=True)
    print(pyclone_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))
    texts = list()

    if args.CCF:
        name = "CCF"
        watching = "cellular_prevalence"
        loc = "lower left"
    elif args.VAF:
        name = "VAF"
        watching = "variant_allele_frequency"
        loc = "upper right"
    else:
        raise Exception("Something went wrong!!")

    matplotlib.pyplot.scatter(pyclone_data.loc[~(pyclone_data["gene_census"]) & ~(pyclone_data["mutation"]) & (pyclone_data["sample_id"] == first_name), watching], pyclone_data.loc[~(pyclone_data["gene_census"]) & ~(pyclone_data["mutation"]) & (pyclone_data["sample_id"] == second_name), watching], c="tab:gray", marker="o", alpha=0.3, s=12 ** 2, edgecolor="none", label="Synonymous mutations")
    matplotlib.pyplot.scatter(pyclone_data.loc[~(pyclone_data["gene_census"]) & (pyclone_data["mutation"]) & (pyclone_data["sample_id"] == first_name), watching], pyclone_data.loc[~(pyclone_data["gene_census"]) & (pyclone_data["mutation"]) & (pyclone_data["sample_id"] == second_name), watching], c="black", marker="*", alpha=0.3, s=15 ** 2, edgecolor="none", label="Functional mutations")
    matplotlib.pyplot.scatter(pyclone_data.loc[(pyclone_data["gene_census"]) & (pyclone_data["mutation"]) & (pyclone_data["sample_id"] == first_name), watching], pyclone_data.loc[(pyclone_data["gene_census"]) & (pyclone_data["mutation"]) & (pyclone_data["sample_id"] == second_name), watching], c="tab:red", marker="*", alpha=1.0, s=20 ** 2, edgecolor="none", label="Cancer genes")

    for mutation in tqdm.tqdm(pyclone_data.loc[(pyclone_data["sample_id"] == first_name) & (pyclone_data["gene_census"]) & (pyclone_data["mutation"]), "mutation_id"]):
        first = pyclone_data.loc[(pyclone_data["mutation_id"] == mutation) & (pyclone_data["sample_id"] == first_name), watching].to_numpy()[0]
        second = pyclone_data.loc[(pyclone_data["mutation_id"] == mutation) & (pyclone_data["sample_id"] == second_name), watching].to_numpy()[0]

        gene = pyclone_data.loc[(pyclone_data["mutation_id"] == mutation) & (pyclone_data["sample_id"] == first_name), "gene"].to_numpy()[0]
        protein = mutation.split(":")[-1]

        if (gene == "nan") or (protein == "nan"):
            continue

        if (first > 0.6) and (second == 0.0):
            texts.append(matplotlib.pyplot.text(first, second, f"{gene}: {protein}", fontsize="xx-small", bbox={"color": "white", "alpha": 0.5}, horizontalalignment="center", verticalalignment="center"))
        elif (first > 0.6) and (second > 0.6):
            texts.append(matplotlib.pyplot.text(first, second, f"{gene}: {protein}", fontsize="xx-small", bbox={"color": "white", "alpha": 0.5}, horizontalalignment="center", verticalalignment="center"))
        elif (first > 0.0) and (second > 0.0):
            texts.append(matplotlib.pyplot.text(first, second, f"{gene}: {protein}", fontsize="xx-small", bbox={"color": "white", "alpha": 0.5}, horizontalalignment="center", verticalalignment="center"))

    matplotlib.pyplot.axline((0, 0), (1, 1), linestyle="--", color="black", alpha=0.3)
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.xlim(-0.1, 1.1)
    matplotlib.pyplot.ylim(-0.1, 1.1)
    matplotlib.pyplot.xlabel(f"{name} of {first_name} ({step00.get_long_sample_type(first_name)})")
    matplotlib.pyplot.ylabel(f"{name} of {second_name} ({step00.get_long_sample_type(second_name)})")
    matplotlib.pyplot.title(f"{first_name} vs. {second_name}")
    matplotlib.pyplot.legend(loc=loc)
    matplotlib.pyplot.tight_layout()
    adjust_text(texts, arrowprops={"arrowstyle": "-", "color": "k", "linewidth": 0.5}, ax=ax, lim=step00.big)

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
