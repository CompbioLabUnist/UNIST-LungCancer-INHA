"""
convert_mutect2_pyclone.py: convert Mutect2 result to PyClone input format
"""
import argparse
import multiprocessing
import tarfile
import pandas
import step00

sequenza_data = pandas.DataFrame()


def find_normal_cn(chromosome: str, start: int, end: int) -> int:
    d = sequenza_data.loc[(sequenza_data["chromosome"] == chromosome) & (sequenza_data["start.pos"] >= start) & (sequenza_data["end.pos"] <= end), "CNt"]

    if d.empty or sum(d.to_numpy()) != 0:
        return 2
    else:
        return 1


def find_minor_cn(chromosome: str, start: int, end: int) -> int:
    d = sequenza_data.loc[(sequenza_data["chromosome"] == chromosome) & (sequenza_data["start.pos"] >= start) & (sequenza_data["end.pos"] <= end), "B"]

    if d.empty or sum(d.to_numpy()) != 0:
        return 1
    else:
        return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("mutect", help="Mutect2 MAF file", type=str)
    parser.add_argument("sequenza", help="Sequenza tar.gz file", type=str)
    parser.add_argument("output", help="Output basename file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if not args.mutect.endswith(".maf"):
        raise ValueError("Mutect must end with .MAF!!")
    elif not args.sequenza.endswith(".tar.gz"):
        raise ValueError("Sequenza must end with .TAR.gz!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    mutect_data = pandas.read_csv(args.mutect, sep="\t", comment="#", low_memory=False)
    mutect_data = mutect_data.loc[(mutect_data["Chromosome"] != "chrX") & (mutect_data["Chromosome"] != "chrY") & (mutect_data["Chromosome"] != "chrM")]
    print(mutect_data)

    with tarfile.open(args.sequenza, "r:gz") as tar:
        txt_file = list(filter(lambda x: x.endswith("_alternative_solutions.txt"), tar.getnames()))[0]
        tar.extract(txt_file, path=step00.tmpfs)
    solution_data = pandas.read_csv(step00.tmpfs + "/" + txt_file, sep="\t")
    with open(args.output + ".Cellularity.txt", "w") as f:
        f.write("%s" % solution_data.loc[0, "cellularity"])

    with tarfile.open(args.sequenza, "r:gz") as tar:
        txt_file = list(filter(lambda x: x.endswith("_segments.txt"), tar.getnames()))[0]
        tar.extract(txt_file, path=step00.tmpfs)
    sequenza_data = pandas.read_csv(step00.tmpfs + "/" + txt_file, sep="\t").dropna(axis="index")
    print(sequenza_data)

    output_data = pandas.DataFrame()
    output_data["mutation_id"] = list(map(lambda x: x[0] + ":" + str(x[1]) + ":" + str(x[2]) + ":" + x[3] + ">" + x[4] + ":" + str(x[5]), zip(mutect_data["Chromosome"], mutect_data["Start_Position"], mutect_data["End_Position"], mutect_data["Reference_Allele"], mutect_data["Tumor_Seq_Allele2"], mutect_data["HGVSp"])))
    output_data["ref_counts"] = list(map(int, mutect_data["t_ref_count"]))
    output_data["var_counts"] = list(map(int, mutect_data["t_alt_count"]))
    with multiprocessing.Pool(args.cpus) as pool:
        output_data["normal_cn"] = pool.starmap(find_normal_cn, zip(mutect_data["Chromosome"], mutect_data["Start_Position"], mutect_data["End_Position"]))
        output_data["minor_cn"] = pool.starmap(find_minor_cn, zip(mutect_data["Chromosome"], mutect_data["Start_Position"], mutect_data["End_Position"]))
    output_data["major_cn"] = output_data["normal_cn"] - output_data["minor_cn"]

    print(output_data)
    output_data.to_csv(args.output + ".VAF.tsv", sep="\t", index=False)
