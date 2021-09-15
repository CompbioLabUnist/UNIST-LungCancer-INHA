"""
convert_mutect2_pyclone.py: convert Mutect2 result to PyClone input format
"""
import argparse
import typing
from adjustText import adjust_text
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import step00

mutenricher_set = set()
census_set = set()


def is_included(gene: str) -> str:
    if (gene in mutenricher_set) and (gene in census_set):
        return "Census+Driver"
    elif gene in mutenricher_set:
        return "Driver"
    elif gene in census_set:
        return "Census"
    else:
        return "None"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("loci", help="Pyclone loci TSV file", type=str)
    parser.add_argument("driver", help="MutEnricher driver gene TSV file (not necessarily TSV)", type=str)
    parser.add_argument("census", help="Cancer gene census CSV file", type=str)
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
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")
    elif not (0 < args.p < 1):
        raise ValueError("P-value must be (0, 1)!!")

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
    print(census_data)

    driver_data = pandas.read_csv(args.driver, sep="\t")
    driver_data = driver_data.loc[(driver_data["Gene"].isin(census_data["Gene Symbol"]))]
    for column in step00.MutEnricher_pval_columns:
        driver_data = driver_data.loc[(driver_data[column] < args.p)]
    driver_data["chromosome"] = list(map(lambda x: x.replace("-", ":").split(":")[0], driver_data["coordinates"]))
    driver_data["start"] = list(map(lambda x: int(x.replace("-", ":").split(":")[1]), driver_data["coordinates"]))
    driver_data["end"] = list(map(lambda x: int(x.replace("-", ":").split(":")[2]), driver_data["coordinates"]))
    print(driver_data)

    gene_names: typing.List[typing.List[str]] = list()
    for index, row in pyclone_data.iterrows():
        tmp_data = driver_data.loc[(driver_data["chromosome"] == row["chromosome"]) & (driver_data["start"] <= row["start"]) & (row["end"] <= driver_data["end"]), "Gene"]

        if tmp_data.empty:
            gene_names.append([])
        else:
            gene_names.append(list(tmp_data.to_numpy()))

    pyclone_data["gene"] = gene_names
    pyclone_data = pyclone_data.explode("gene", ignore_index=True)
    pyclone_data.dropna(inplace=True)
    print(pyclone_data)

    mutenricher_set = set(driver_data["Gene"])
    print(driver_data)

    pyclone_data["Database"] = list(map(is_included, pyclone_data["gene"]))
    pyclone_data.rename(columns={"cluster_id": "Cluster ID"}, inplace=True)
    print(pyclone_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
    texts = list()
    if args.CCF:
        seaborn.lineplot(data=pyclone_data, x="sample_id", y="cellular_prevalence", hue="Cluster ID", style="Database", legend="brief", ax=ax, estimator=None, units="gene")
        matplotlib.pyplot.ylabel("Cancer Cell Fraction")
        for index, row in pyclone_data.iterrows():
            if row["Database"] == "Census+Driver":
                texts.append(matplotlib.pyplot.text(row["sample_id"], row["cellular_prevalence"], row["gene"], fontsize="small", horizontalalignment="center", bbox={"facecolor": "white", "alpha": 0.5}))
    elif args.VAF:
        seaborn.lineplot(data=pyclone_data, x="sample_id", y="variant_allele_frequency", hue="Cluster ID", style="Database", legend="brief", ax=ax, estimator=None, units="gene")
        matplotlib.pyplot.ylabel("Variant Allele Frequency")
        for index, row in pyclone_data.iterrows():
            if row["Database"] == "Census+Driver":
                texts.append(matplotlib.pyplot.text(row["sample_id"], row["variant_allele_frequency"], row["gene"], fontsize="small", horizontalalignment="center", bbox={"facecolor": "white", "alpha": 0.5}))
    else:
        raise Exception("Something went wrong!!")

    matplotlib.pyplot.xlabel("Samples")
    matplotlib.pyplot.ylim(-0.1, 1.1)
    adjust_text(texts, arrowprops={"arrowstyle": "-", "color": "k"}, ax=ax, lim=10 ** 6)

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
