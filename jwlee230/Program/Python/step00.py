"""
step00.py: for base implementation
"""
import os
import re
import typing
import numpy
import pandas

tmpfs = "/tmpfs"
chromosome_list = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22".split()
chromosome_full_list = chromosome_list + ["chrX"]
big = 10 ** 6
small = 10 ** 3
WES_length = 36995809

matplotlib_parameters = {"font.size": 50, "axes.labelsize": 50, "axes.titlesize": 75, "xtick.labelsize": 50, "ytick.labelsize": 50, "font.family": "Arial", "legend.fontsize": 30, "legend.title_fontsize": 30, "figure.dpi": 200, "pdf.fonttype": 42, "ps.fonttype": 42}

long_sample_type_dict = {"N": "Normal", "C": "CIS+AIS", "A": "AAH", "P": "Primary", "D": "Dysplasia", "M": "MIA"}
long_sample_type_list = ["Normal", "Dysplasia", "AAH", "CIS+AIS", "MIA", "Primary"]
sample_order_dict = {"Normal": 0, "Dysplasia": 1, "CIS+AIS": 2, "AAH": 1, "MIA": 3, "Primary": 4}
simple_stage_list = ("Normal", "Precancer", "Primary")

stage_color_code = {"Normal": "cyan", "Dysplasia": "olive", "AAH": "orange", "CIS+AIS": "red", "MIA": "brown", "Primary": "silver"}
stage_linestyle = {"Normal": "-", "Dysplasia": ":", "AAH": ":", "CIS+AIS": "-.", "MIA": "--", "Primary": "-"}

MutEnricher_pval_columns = ["Fisher_FDR", "Fisher_pval", "gene_pval"]

nonsynonymous_column = "Variant_Classification"
nonsynonymous_mutations = {"Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation"}
nonsynonymous_notations = {"Nonsense_Mutation": "Nonsense", "In_Frame_Del": "In frame indel", "In_Frame_Ins": "In frame indel", "Frame_Shift_Del": "Frameshift indel", "Missense_Mutation": "Missense", "Splice_Site": "Splice site", "Frame_Shift_Ins": "Frameshift indel", "Translation_Start_Site": "TSS", "Nonstop_Mutation": "Nonstop"}
nonsynonymous_coloring = {"Missense": "darkgreen", "Nonsense": "cyan", "In frame indel": "navy", "Frameshift indel": "gold", "Splice site": "darkviolet", "LOH": "orange", "TSS": "chocolate", "Nonstop": "crimson", "Absent": "lightgray", "Multiple": "black"}

cnv_mutations = nonsynonymous_mutations | {"CNV loss", "CNV gain"}
cnv_coloring = nonsynonymous_coloring | {"CNV loss": {"facecolor": "none", "edgecolor": "tab:cyan", "linewidth": 9}, "CNV gain": {"facecolor": "none", "edgecolor": "tab:red", "linewidth": 9}}

derivations = ("Accuracy", "Balanced Accuracy", "Sensitivity", "Specificity", "Precision")
venn_format = "{size:d} ({percentage:.1f}%)"
Gistic_columns = ["Sample", "Chromosome", "Start Position", "End Position", "Num markers", "Seg.CN"]

Tcell_list = ["CD4+ Th", "CD8 low T", "CD8+/CD4+ Mixed Th", "Cytotoxic CD8+ T", "Exhausted CD8+ T", "Exhausted Tfh", "Naive CD4+ T", "Naive CD8+ T", "Treg"]

PathSeq_type_list = ["class", "family", "genus", "kingdom", "no_rank", "order", "phylum", "root", "species", "species_group", "species_subgroup", "subclass", "subfamily", "subgenus", "subkingdom", "suborder", "subphylum", "subspecies", "superkingdom", "tribe", "varietas"]

sharing_strategy = ["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]
sharing_columns = ["Mutation Shared Proportion", "Mutation Shared Proportion (SYN)", "Mutation Shared Proportion (Union)", "Mutation Shared Proportion (Union & SYN)", "Mutation Shared Count", "Mutation Shared Count (SYN)", "Mutation Shared Count per TMB", "Mutation Shared Count (SYN) per TMB"]


def file_list(path: str) -> typing.List[str]:
    """
    file_list: return a list of files in path
    """
    return list(filter(lambda x: os.path.isfile(x), list(map(lambda x: os.path.join(path, x), os.listdir(path)))))


def directory_list(path: str) -> typing.List[str]:
    """
    directory_list: return a list of directories in path
    """
    return list(filter(lambda x: os.path.isdir(x), list(map(lambda x: os.path.join(path, x), os.listdir(path)))))


def get_id(ID: str) -> str:
    return ID.split("/")[-1].split(".")[0]


def get_patient(ID: str) -> str:
    """
    get_patient: get patient ID from sample ID
    """
    return re.findall(r"(^(cn)?\d+)", get_id(ID))[0][0]


def get_paired_normal(ID: str) -> str:
    """
    get_paired_normal: get paired normal sample ID
    """
    return get_patient(ID) + "N"


def get_paired_primary(ID: str) -> str:
    ID = get_id(ID)
    if ID in ["cn107A1"]:
        return ID[:-2] + "P"
    elif ID in ["12N", "14N"]:
        return ID[:-1] + "P1"
    elif ID[-1].isalpha():
        return ID[:-1] + "P"
    else:
        return ID[:-2] + "P1"


def get_sample_type(ID: str) -> str:
    """
    get_sample_type: get the type of sample from sample ID
    """
    return re.findall(r"[PCADMN]", get_id(ID))[0][0]


def get_long_sample_type(ID: str) -> str:
    """
    get_long_sample_type: get the full type of sample from sample ID
    """
    return long_sample_type_dict[get_sample_type(ID)]


def get_simple_sample_type(ID: str) -> str:
    """
    get_simple_sample_type: similar to get_sample_type, but only pre-cancer
    """
    t = get_long_sample_type(ID)
    if t in ["Normal", "Primary"]:
        return t
    else:
        return "Precancer"


def list_first_last(li: typing.List[typing.Any], el: typing.Any) -> typing.Tuple[int, int]:
    """
    list_first_last: get the first and last one which specified
    """
    return li.index(el), max(loc for loc, val in enumerate(li) if val == el)


def sorting(ID: str) -> typing.Tuple[str, int, str]:
    """
    sorting: sorting key by patient-type
    """
    ID = get_id(ID)
    return (get_patient(ID), long_sample_type_list.index(get_long_sample_type(ID)), ID)


def sorting_by_type(ID: str) -> typing.Tuple[int, str, str]:
    """
    sorting_by_type: sorting key by type-sample
    """
    ID = get_id(ID)
    return (long_sample_type_list.index(get_long_sample_type(ID)), get_patient(ID), ID)


def get_color_by_type(ID: str) -> str:
    """
    get_color_by_type: get color by type
    """
    return stage_color_code[get_long_sample_type(ID)]


def get_clinical_data(filename: str) -> pandas.DataFrame:
    """
    get_clinical_data: get clinical data for select proper patients
    """
    return pandas.read_csv(filename, index_col="Serial_No", skiprows=[1], verbose=True).dropna(axis="index", how="all")


def aggregate_confusion_matrix(confusion_matrix: numpy.ndarray, derivation: str = "") -> float:
    """
    aggregate_confusion_matrix: derivations from confusion matrix
    """

    assert (derivation in derivations)
    assert confusion_matrix.shape == (2, 2)

    TP, FP, FN, TN = confusion_matrix[0][0], confusion_matrix[0][1], confusion_matrix[1][0], confusion_matrix[1][1]
    assert TP and FP and FN and TN

    if derivation == "Sensitivity":
        return TP / (TP + FN)
    elif derivation == "Specificity":
        return TN / (TN + FP)
    elif derivation == "Precision":
        return TP / (TP + FP)
    elif derivation == "Accuracy":
        return (TP + TN) / (TP + TN + FP + FN)
    elif derivation == "Balanced Accuracy":
        return TP / (2 * (TP + FN)) + TN / (2 * (TN + FP))
    else:
        raise Exception("Something went wrong!!")


def get_band_data(filename: str) -> pandas.DataFrame:
    """
    get_band_data: get chromosome band data including p-arm and q-arm information
    """
    data = pandas.read_csv(filename, sep="\t", names=["chrom", "chrom_start", "chrom_end", "name", "gie_stain"], verbose=True).dropna(axis="index")
    data["arm"] = list(map(lambda x: x[0], data["name"]))
    return data
