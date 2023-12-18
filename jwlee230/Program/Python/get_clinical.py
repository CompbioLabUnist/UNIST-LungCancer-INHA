"""
get_clinical.py: get clinical data with tableon package
"""
import argparse
import tarfile
import numpy
import pandas
import tableone
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Mutect2 MAF file(s)", type=str, nargs="+")
    parser.add_argument("clinical", help="Clinidata data TSV file w/ mutation shared proportion", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--percentage", help="Percentage of genes to draw", type=float, default=0.1)

    group_subtype = parser.add_mutually_exclusive_group(required=True)
    group_subtype.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group_subtype.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".maf"), args.input)):
        raise ValueError("Input must end with .MAF!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif not (0.0 < args.percentage < 0.5):
        raise ValueError("Percentage must be between 0.0 and 0.5!!")

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(len(patients), sorted(patients))

    patients &= set(map(step00.get_patient, list(filter(lambda x: step00.get_long_sample_type(x) == "Primary", args.input))))
    print(len(patients), sorted(patients))

    clinical_data = clinical_data.loc[sorted(patients), :]
    print(clinical_data)

    categorical = ["Gender", "Smoking-Detail", "pT", "pN", "Stage"]
    continuous = ["Pack-Year", "Overall Survival", "Recurrence-Free Survival"]

    files = list()
    table = tableone.TableOne(clinical_data, columns=categorical + continuous, categorical=categorical)

    files.append("clinical.tsv")
    table.to_csv(files[-1], sep="\t")

    files.append("clinical.xlsx")
    table.to_excel(files[-1])

    for MSP in tqdm.tqdm(step00.sharing_columns):
        lower, higher = numpy.quantile(clinical_data[MSP], args.percentage), numpy.quantile(clinical_data[MSP], 1 - args.percentage)

        drawing_data = clinical_data.copy()
        drawing_data[MSP] = list(map(lambda x: ("Lower") if (x <= lower) else (("Higher") if (x >= higher) else None), clinical_data[MSP]))
        drawing_data = drawing_data.dropna(axis="index", subset=[MSP])

        table = tableone.TableOne(drawing_data, columns=categorical + continuous, categorical=categorical, groupby=MSP, pval=True, htest_name=True)

        files.append(f"{MSP}-percentage.tsv")
        table.to_csv(files[-1], sep="\t")

        files.append(f"{MSP}-percentage.xlsx")
        table.to_excel(files[-1])

    for MSP in tqdm.tqdm(step00.sharing_columns):
        threshold = numpy.median(clinical_data[MSP])

        drawing_data = clinical_data.copy()
        drawing_data[MSP] = list(map(lambda x: "Lower" if (x <= threshold) else "Higher", clinical_data[MSP]))

        table = tableone.TableOne(drawing_data, columns=categorical + continuous, categorical=categorical, groupby=MSP, pval=True, htest_name=True)

        files.append(f"{MSP}-median.tsv")
        table.to_csv(files[-1], sep="\t")

        files.append(f"{MSP}-median.xlsx")
        table.to_excel(files[-1])

    with tarfile.open(args.output, "w") as tar:
        for file in tqdm.tqdm(files):
            tar.add(file, arcname=file)
