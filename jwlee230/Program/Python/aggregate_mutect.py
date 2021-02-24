"""
aggregate_mutect.py: aggregate mutect .MAF files
"""
import argparse
import re
import pandas
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 input .MAF files", type=str, nargs="+")
    parser.add_argument("output", help="Output file", type=str)
    parser.add_argument("ID", help="Patient ID", type=str)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("INPUT must end with .MAF!!")

    args.input.sort()
    IDs = list(map(lambda x: x.split(".")[0], list(map(lambda x: x.split("/")[-1], args.input))))

    raw_data = list()
    for ID, file_name in zip(IDs, args.input):
        if re.findall(r"(^(cn)?\d+)", ID)[0][0] != args.ID:
            continue

        print(ID)
        data = pandas.read_csv(file_name, sep="\t", comment="#", dtype=str)
        data["ID"] = ID
        data["Patient"] = re.findall(r"(^(cn)?\d+)", ID)[0][0]
        raw_data.append(data)

    mutect_data = pandas.concat(raw_data, ignore_index=True, verify_integrity=True)
    del raw_data
    print(mutect_data)
    mutect_data.info()

    step00.make_pickle(args.output, mutect_data)
