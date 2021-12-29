"""
get_gene_location.py: get Genome Location by chromosome
"""
import argparse
import itertools
import multiprocessing
import tarfile
import typing
import pandas
import tqdm
import step00

band_data = pandas.DataFrame()


def get_chromosome(location: str) -> str:
    return "chr" + location.split(":")[0]


def get_start(location: str) -> int:
    return int(location.replace("-", ":").split(":")[1])


def get_end(location: str) -> int:
    return int(location.replace("-", ":").split(":")[2])


def query_band(chromosome: str, start: int, end: int) -> typing.List[str]:
    answer = list()

    answer += list(band_data.loc[(band_data["chrom"] == chromosome) & (band_data["chrom_start"] <= start) & (start <= band_data["chrom_end"]), "chrom-arm"])
    answer += list(band_data.loc[(band_data["chrom"] == chromosome) & (start <= band_data["chrom_start"]) & (band_data["chrom_end"] <= end), "chrom-arm"])
    answer += list(band_data.loc[(band_data["chrom"] == chromosome) & (band_data["chrom_start"] <= end) & (end <= band_data["chrom_end"]), "chrom-arm"])

    return sorted(set(answer))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("one", help="CGC Tier 1 CSV files", type=str)
    parser.add_argument("both", help="CGC Both tier CSV files", type=str)
    parser.add_argument("band", help="Chromosome band txt file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if not args.one.endswith(".csv"):
        raise ValueError("One must end with .CSV!!")
    elif not args.both.endswith(".csv"):
        raise ValueError("Both must end with .CSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    chromosome_list = list(map(lambda x: "-".join(x), itertools.product(step00.chromosome_list, ["p", "q"])))
    file_list = list()
    length_limit = 5

    band_data = step00.get_band_data(args.band)
    band_data["chrom-arm"] = list(map(lambda x: "-".join(x), band_data[["chrom", "arm"]].itertuples(index=False, name=None)))
    print(band_data)

    cgc_one_data = pandas.read_csv(args.one)
    cgc_one_data = cgc_one_data.loc[~(cgc_one_data["Genome Location"].str.contains(":-"))]
    with multiprocessing.Pool(args.cpus) as pool:
        cgc_one_data["Chromosome"] = pool.map(get_chromosome, cgc_one_data["Genome Location"])
        cgc_one_data["Start"] = pool.map(get_start, cgc_one_data["Genome Location"])
        cgc_one_data["End"] = pool.map(get_end, cgc_one_data["Genome Location"])
        cgc_one_data["Chromosome-Arm"] = pool.starmap(query_band, cgc_one_data[["Chromosome", "Start", "End"]].itertuples(index=False, name=None))
    cgc_one_data = cgc_one_data.explode("Chromosome-Arm", ignore_index=True)
    cgc_one_set = set(cgc_one_data["Gene Symbol"])
    print(cgc_one_data)

    cgc_both_data = pandas.read_csv(args.both)
    cgc_both_data = cgc_both_data.loc[~(cgc_both_data["Genome Location"].str.contains(":-"))]
    with multiprocessing.Pool(args.cpus) as pool:
        cgc_both_data["Chromosome"] = pool.map(get_chromosome, cgc_both_data["Genome Location"])
        cgc_both_data["Start"] = pool.map(get_start, cgc_both_data["Genome Location"])
        cgc_both_data["End"] = pool.map(get_end, cgc_both_data["Genome Location"])
        cgc_both_data["Chromosome-Arm"] = pool.starmap(query_band, cgc_both_data[["Chromosome", "Start", "End"]].itertuples(index=False, name=None))
    cgc_both_data = cgc_both_data.explode("Chromosome-Arm", ignore_index=True)
    print(cgc_both_data)

    for chromosome in tqdm.tqdm(chromosome_list):
        tmp_data = cgc_one_data.loc[(cgc_one_data["Chromosome-Arm"] == chromosome), ["Gene Symbol", "Name", "Chromosome", "Start", "End"]]
        if tmp_data.empty:
            tmp_data = pandas.DataFrame(data=[["None", "", "", "", ""]], columns=["Gene Symbol", "Name", "Chromosome", "Start", "End"])

        file_list.append("{0}.CGC-1.tsv".format(chromosome))
        tmp_data.to_csv(file_list[-1], sep="\t", header=True, index=False)

        tmp_data = tmp_data.loc[:, ["Gene Symbol", "Name"]]
        if (l := tmp_data.shape[0]) > length_limit:
            tmp_data.columns = ["Gene Symbol ({0})".format(l), "Name"]
            tmp_data = tmp_data.iloc[:length_limit, :]

        file_list.append("{0}.CGC-1.tex".format(chromosome))
        tmp_data.to_latex(file_list[-1], header=True, index=False)

    for chromosome in tqdm.tqdm(chromosome_list):
        tmp_data = cgc_both_data.loc[(cgc_both_data["Chromosome-Arm"] == chromosome) & ~(cgc_both_data["Gene Symbol"].isin(cgc_one_set)), ["Gene Symbol", "Name", "Chromosome", "Start", "End"]]
        if tmp_data.empty:
            tmp_data = pandas.DataFrame(data=[["None", "", "", "", ""]], columns=["Gene Symbol", "Name", "Chromosome", "Start", "End"])

        file_list.append("{0}.CGC-2.tsv".format(chromosome))
        tmp_data.to_csv(file_list[-1], sep="\t", header=True, index=False)

        tmp_data = tmp_data.loc[:, ["Gene Symbol", "Name"]]
        if (l := tmp_data.shape[0]) > length_limit:
            tmp_data.columns = ["Gene Symbol ({0})".format(l), "Name"]
            tmp_data = tmp_data.iloc[:length_limit, :]

        file_list.append("{0}.CGC-2.tex".format(chromosome))
        tmp_data.to_latex(file_list[-1], header=True, index=False)

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(file_list):
            tar.add(f, arcname=f)
