"""
count_SharedMutation.py: count Shared Mutation
"""
import argparse
import multiprocessing
import pandas
import tqdm
import step00


def read_maf(filename: str) -> pandas.DataFrame:
    return pandas.read_csv(filename, sep="\t", comment="#", low_memory=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF file(s)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinidata data CSV file", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")
    elif not args.clinical.endswith(".csv"):
        raise ValueError("Clinical must end with .CSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    clinical_data = step00.get_clinical_data(args.clinical)
    print(clinical_data)

    patients = set(clinical_data.index)
    print("Patients:", len(patients))

    args.input = list(filter(lambda x: step00.get_patient(x) in patients, args.input))
    args.input.sort(key=step00.sorting)
    with multiprocessing.Pool(args.cpus) as pool:
        mutect_data = pandas.concat(pool.map(read_maf, args.input), ignore_index=True, copy=False)
        mutect_data["Tumor_Sample_Barcode"] = pool.map(step00.get_id, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Patient"] = pool.map(step00.get_patient, mutect_data["Tumor_Sample_Barcode"])
        mutect_data["Stage"] = pool.map(step00.get_long_sample_type, mutect_data["Tumor_Sample_Barcode"])
    print(mutect_data)

    precancer_samples = list(filter(lambda x: step00.get_long_sample_type(x) not in {"Normal", "Primary"}, mutect_data["Tumor_Sample_Barcode"].unique()))

    raw_output_data = list()
    for precancer_sample in tqdm.tqdm(precancer_samples):
        primary_sample = step00.get_paired_primary(precancer_sample)
        patient = step00.get_patient(precancer_sample)

        precancer_set = set(mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == precancer_sample), ["Hugo_Symbol"] + step00.sharing_strategy + [step00.nonsynonymous_column]].itertuples(index=False, name=None))
        primary_set = set(mutect_data.loc[(mutect_data["Tumor_Sample_Barcode"] == primary_sample), ["Hugo_Symbol"] + step00.sharing_strategy + [step00.nonsynonymous_column]].itertuples(index=False, name=None))

        raw_output_data += list(map(lambda x: x + (precancer_sample, primary_sample, patient), list(precancer_set & primary_set)))

    output_data = pandas.DataFrame(raw_output_data, columns=["Hugo_Symbol"] + step00.sharing_strategy + [step00.nonsynonymous_column, "Precancer", "Primary", "Patient"]).sort_values(["Start_Position", "End_Position"]).sort_values("Chromosome", key=lambda x: x.apply(lambda y: (step00.chromosome_list.index(y), y) if (y in step00.chromosome_list) else (len(step00.chromosome_list), y)), ignore_index=True, kind="mergesort")
    output_data.index.name = "Index"
    print(output_data)
    output_data.to_csv(args.output, sep="\t", index=True)
