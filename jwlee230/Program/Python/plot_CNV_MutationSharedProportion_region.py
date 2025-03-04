"""
plot_CNV_MutationSharedProportion_region.py: violin plot & joint plot Copy Number Variation with Mutation Shared Proportion for region length
"""
import argparse
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import statannotations.Annotator
import tqdm
import step00

input_data = pandas.DataFrame()
watching = ""
threshold = 0.0


def get_chromosome_data(sample: str) -> int:
    tmp_data = input_data.loc[(input_data["Sample"] == sample) & (((input_data[watching] * input_data["weight"]) <= (1 - threshold)) | ((input_data[watching] * input_data["weight"]) >= (1 + args.threshold))), :]
    return (sum(tmp_data["end"]) - sum(tmp_data["start"]) + tmp_data.shape[0]) / (10 ** 9)


def get_chromosome_data_loss(sample: str) -> int:
    tmp_data = input_data.loc[(input_data["Sample"] == sample) & ((input_data[watching] * input_data["weight"]) <= (1 - threshold)), :]
    return (sum(tmp_data["end"]) - sum(tmp_data["start"]) + tmp_data.shape[0]) / (10 ** 9)


def get_chromosome_data_gain(sample: str) -> int:
    tmp_data = input_data.loc[(input_data["Sample"] == sample) & ((input_data[watching] * input_data["weight"]) >= (1 + threshold)), :]
    return (sum(tmp_data["end"]) - sum(tmp_data["start"]) + tmp_data.shape[0]) / (10 ** 9)


def get_middle(values):
    return (min(values) + max(values)) / 2


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="CNV segment.tsv file", type=str)
    parser.add_argument("clinical", help="Clinical data with Mutation Shared Proportion TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--watching", help="Watching column name", type=str, required=True)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)
    parser.add_argument("--threshold", help="Threshold for gain/loss", type=float, default=0.2)
    parser.add_argument("--percentage", help="Percentage of patients to include", type=float, default=0.1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--SQC", help="Get SQC patient only", action="store_true", default=False)
    group.add_argument("--ADC", help="Get ADC patient only", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.clinical.endswith(".tsv"):
        raise ValueError("Clinical must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")
    elif not (0 < args.threshold < 1):
        raise ValueError("Threshold must be (0, 1)")
    elif not (0.0 < args.percentage < 0.5):
        raise ValueError("Percentage must be (0.0, 0.5)!!")

    watching = args.watching
    threshold = args.threshold

    clinical_data = pandas.read_csv(args.clinical, sep="\t", index_col=0)
    print(clinical_data)

    if args.SQC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "SQC")].index)
    elif args.ADC:
        patients = set(clinical_data.loc[(clinical_data["Histology"] == "ADC")].index)
    else:
        raise Exception("Something went wrong!!")
    print(len(patients), patients)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data = input_data.loc[(input_data["Patient"].isin(patients))]
    print(input_data)

    clinical_data = clinical_data.loc[list(filter(lambda x: x in patients, list(clinical_data.index)))]
    print(clinical_data)

    sample_list = sorted(set(input_data["Sample"]), key=step00.sorting_by_type)
    print(sample_list)

    output_data = pandas.DataFrame(data=sample_list, columns=["Sample"])
    with multiprocessing.Pool(args.cpus) as pool:
        output_data["Patient"] = pool.map(step00.get_patient, output_data["Sample"])
        output_data["Stage"] = pool.map(step00.get_long_sample_type, output_data["Sample"])
        output_data["Region"] = pool.map(get_chromosome_data, output_data["Sample"])
        output_data["Region-Loss"] = pool.map(get_chromosome_data_loss, output_data["Sample"])
        output_data["Region-Gain"] = pool.map(get_chromosome_data_gain, output_data["Sample"])
    print(output_data)

    sample_list = sorted(set(output_data["Sample"]), key=step00.sorting_by_type)
    print(sample_list)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    MSP_order = ["PSM-L", "PSM-H"]
    compare_list = [(("PSM-L", "Precancer"), ("PSM-L", "Primary")), (("PSM-H", "Precancer"), ("PSM-H", "Primary")), (("PSM-L", "Precancer"), ("PSM-H", "Precancer")), (("PSM-L", "Primary"), ("PSM-H", "Primary"))]

    for MSP in tqdm.tqdm(step00.sharing_columns[:2]):
        lower_bound, higher_bound = numpy.quantile(clinical_data[MSP], args.percentage), numpy.quantile(clinical_data[MSP], 1 - args.percentage)

        drawing_data = output_data.loc[(output_data["Sample"].isin(clinical_data[f"{MSP}-sample"])) | (output_data["Sample"].isin(list(map(step00.get_paired_primary, clinical_data[f"{MSP}-sample"]))))].copy()
        drawing_data[MSP] = list(map(lambda x: "PSM-L" if (clinical_data.loc[x, MSP] <= lower_bound) else ("PSM-H" if (clinical_data.loc[x, MSP] >= higher_bound) else "NS"), drawing_data["Patient"]))
        drawing_data["Stage"] = list(map(lambda x: "Primary" if (x == "Primary") else "Precancer", drawing_data["Stage"]))
        drawing_data = drawing_data.loc[(drawing_data[MSP].isin(MSP_order))]

        stage_list = ["Precancer", "Primary"]

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.violinplot(data=drawing_data, x=MSP, order=MSP_order, y="Region", hue="Stage", hue_order=stage_list, palette=step00.precancer_color_code, inner="box", cut=0, linewidth=5, ax=ax)
        statannotations.Annotator.Annotator(ax, compare_list, data=drawing_data, x=MSP, order=MSP_order, y="Region", hue="Stage", hue_order=stage_list).configure(test="Mann-Whitney", text_format="star", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()

        matplotlib.pyplot.xlabel("")
        matplotlib.pyplot.ylabel("Size of somatic CNV region (Gbp)")
        matplotlib.pyplot.tight_layout()

        figures.append(f"Violin_{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

        fig, axs = matplotlib.pyplot.subplots(figsize=(18 * 2, 18), ncols=2)

        seaborn.violinplot(data=drawing_data, x=MSP, order=MSP_order, y="Region-Loss", hue="Stage", hue_order=stage_list, palette=step00.precancer_color_code, inner="box", cut=0, linewidth=5, ax=axs[0])
        statannotations.Annotator.Annotator(axs[0], compare_list, data=drawing_data, x=MSP, order=MSP_order, y="Region-Loss", hue="Stage", hue_order=stage_list).configure(test="Mann-Whitney", text_format="star", loc="inside", verbose=0).apply_and_annotate()

        seaborn.violinplot(data=drawing_data, x=MSP, order=MSP_order, y="Region-Gain", hue="Stage", hue_order=stage_list, palette=step00.precancer_color_code, inner="box", cut=0, ax=axs[1])
        statannotations.Annotator.Annotator(axs[1], compare_list, data=drawing_data, x=MSP, order=MSP_order, y="Region-Gain", hue="Stage", hue_order=stage_list).configure(test="Mann-Whitney", text_format="star", loc="inside", verbose=0, comparisons_correction=None).apply_and_annotate()

        axs[0].set_title("CNV Loss")
        axs[1].set_title("CNV Gain")
        axs[0].set_xlabel("")
        axs[1].set_xlabel("")
        axs[0].set_ylabel("Size of somatic CNV region (Gbp)")
        axs[1].set_ylabel("Size of somatic CNV region (Gbp)")
        matplotlib.pyplot.tight_layout()

        figures.append(f"Violin_LossGain_{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    stage_list = list(filter(lambda x: output_data.loc[(output_data["Stage"] == x)].shape[0] > 3, step00.long_sample_type_list))

    output_data = output_data.loc[(output_data["Stage"].isin(stage_list))]

    for MSP in tqdm.tqdm(step00.sharing_columns[:2]):
        output_data[MSP] = list(map(lambda x: clinical_data.loc[x, MSP], output_data["Patient"]))

        r, p = scipy.stats.spearmanr(output_data[MSP], output_data["Region"])
        g = seaborn.jointplot(data=output_data, x=MSP, y="Region", hue="Stage", hue_order=stage_list, palette=step00.stage_color_code, height=18, ratio=5, kind="scatter")
        if p < 1e-3:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        g.set_axis_labels(MSP, "Number of somatic CNV region (Gbp)")
        figures.append(f"Joint_All_{MSP}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))
        seaborn.regplot(data=output_data, x=MSP, y="Region", fit_reg=True, scatter=True, truncate=False, ax=ax)

        if p < 1e-3:
            matplotlib.pyplot.text(get_middle(output_data[MSP]), get_middle(output_data["Region"]), f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            matplotlib.pyplot.text(get_middle(output_data[MSP]), get_middle(output_data["Ploidy"]), f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        matplotlib.pyplot.xlabel("PSM")
        matplotlib.pyplot.ylabel("Number of somatic CNV region (Gbp)")
        matplotlib.pyplot.tight_layout()

        figures.append(f"Scatter_All_{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

        r, p = scipy.stats.spearmanr(output_data[MSP], output_data["Region-Loss"])
        g = seaborn.jointplot(data=output_data, x=MSP, y="Region-Loss", hue="Stage", hue_order=stage_list, palette=step00.stage_color_code, height=18, ratio=5, kind="scatter")
        if p < 1e-3:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        g.set_axis_labels(MSP, "Number of somatic CNV-Loss region (Gbp)")
        figures.append(f"Joint_All-Loss_{MSP}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.regplot(data=output_data, x=MSP, y="Region-Loss", fit_reg=True, scatter=True, truncate=False, ax=ax)

        if p < 1e-3:
            matplotlib.pyplot.text(get_middle(output_data[MSP]), get_middle(output_data["Region-Loss"]), f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            matplotlib.pyplot.text(get_middle(output_data[MSP]), get_middle(output_data["Region-Loss"]), f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.xlabel("PSM")
        matplotlib.pyplot.ylabel("Number of somatic CNV-Loss region (Gbp)")
        matplotlib.pyplot.tight_layout()

        figures.append(f"Scatter_All-Loss_{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

        r, p = scipy.stats.spearmanr(output_data[MSP], output_data["Region-Gain"])
        g = seaborn.jointplot(data=output_data, x=MSP, y="Region-Gain", hue="Stage", hue_order=stage_list, palette=step00.stage_color_code, height=18, ratio=5, kind="scatter")
        if p < 1e-3:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        g.set_axis_labels(MSP, "Number of somatic CNV-Gain region (Gbp)")
        figures.append(f"Joint_All-Gain_{MSP}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.regplot(data=output_data, x=MSP, y="Region-Gain", fit_reg=True, scatter=True, truncate=False, ax=ax)

        if p < 1e-3:
            matplotlib.pyplot.text(get_middle(output_data[MSP]), get_middle(output_data["Region-Gain"]), f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            matplotlib.pyplot.text(get_middle(output_data[MSP]), get_middle(output_data["Region-Gain"]), f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.xlabel("PSM")
        matplotlib.pyplot.ylabel("Number of somatic CNV-Gain region (Gbp)")
        matplotlib.pyplot.tight_layout()

        figures.append(f"Scatter_All-Gain_{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    for stage, MSP in tqdm.tqdm(list(itertools.product(stage_list, step00.sharing_columns[:2]))):
        tmp_data = output_data.loc[(output_data["Stage"] == stage)]

        r, p = scipy.stats.spearmanr(tmp_data[MSP], tmp_data["Region"])
        g = seaborn.jointplot(data=tmp_data, x=MSP, y="Region", color=step00.stage_color_code[stage], height=18, ratio=5, kind="reg")
        if p < 1e-3:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        g.set_axis_labels(MSP, "Number of somatic CNV region (Gbp)")
        figures.append(f"Joint_{stage}_{MSP}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.regplot(data=tmp_data, x=MSP, y="Region", color=step00.stage_color_code[stage], fit_reg=True, scatter=True, truncate=False, ax=ax)

        if p < 1e-3:
            matplotlib.pyplot.text(get_middle(tmp_data[MSP]), get_middle(tmp_data["Region"]), f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            matplotlib.pyplot.text(get_middle(tmp_data[MSP]), get_middle(tmp_data["Region"]), f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.xlabel("PSM")
        matplotlib.pyplot.ylabel("Number of somatic CNV region (Gbp)")
        matplotlib.pyplot.tight_layout()

        figures.append(f"Scatter_{stage}_{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

        r, p = scipy.stats.spearmanr(tmp_data[MSP], tmp_data["Region-Loss"])
        g = seaborn.jointplot(data=tmp_data, x=MSP, y="Region-Loss", color=step00.stage_color_code[stage], height=18, ratio=5, kind="reg")
        if p < 1e-3:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        g.set_axis_labels(MSP, "Number of somatic CNV-Loss region (Gbp)")
        figures.append(f"Joint_{stage}-Loss_{MSP}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.regplot(data=tmp_data, x=MSP, y="Region-Loss", color=step00.stage_color_code[stage], fit_reg=True, scatter=True, truncate=False, ax=ax)
        if p < 1e-3:
            matplotlib.pyplot.text(get_middle(tmp_data[MSP]), get_middle(tmp_data["Region-Loss"]), f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            matplotlib.pyplot.text(get_middle(tmp_data[MSP]), get_middle(tmp_data["Region-Loss"]), f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.xlabel("PSM")
        matplotlib.pyplot.ylabel("Number of somatic CNV-Loss region (Gbp)")
        matplotlib.pyplot.tight_layout()

        figures.append(f"Scatter_{stage}-Loss_{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

        r, p = scipy.stats.spearmanr(tmp_data[MSP], tmp_data["Region-Gain"])
        g = seaborn.jointplot(data=tmp_data, x=MSP, y="Region-Gain", color=step00.stage_color_code[stage], height=18, ratio=5, kind="reg")
        if p < 1e-3:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        g.set_axis_labels(MSP, "Number of somatic CNV-Gain region (Gbp)")
        figures.append(f"Joint_{stage}-Gain_{MSP}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.regplot(data=tmp_data, x=MSP, y="Region-Gain", color=step00.stage_color_code[stage], fit_reg=True, scatter=True, truncate=False, ax=ax)

        if p < 1e-3:
            matplotlib.pyplot.text(get_middle(tmp_data[MSP]), get_middle(tmp_data["Region-Gain"]), f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            matplotlib.pyplot.text(get_middle(tmp_data[MSP]), get_middle(tmp_data["Region-Gain"]), f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.xlabel("PSM")
        matplotlib.pyplot.ylabel("Number of somatic CNV-Gain region (Gbp)")
        matplotlib.pyplot.tight_layout()

        figures.append(f"Scatter_{stage}-Gain_{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    stage_list = ["Precancer", "Primary"]

    output_data["Stage"] = list(map(lambda x: "Primary" if (x == "Primary") else "Precancer", output_data["Stage"]))
    for MSP in tqdm.tqdm(step00.sharing_columns[:2]):
        output_data[MSP] = list(map(lambda x: clinical_data.loc[x, MSP], output_data["Patient"]))

        r, p = scipy.stats.spearmanr(output_data[MSP], output_data["Region"])
        g = seaborn.jointplot(data=output_data, x=MSP, y="Region", hue="Stage", hue_order=stage_list, palette=step00.precancer_color_code, height=18, ratio=5, kind="scatter")
        if p < 1e-3:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        g.set_axis_labels(MSP, "Number of somatic CNV region (Gbp)")
        figures.append(f"Joint_Precancer_{MSP}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

        r, p = scipy.stats.spearmanr(output_data[MSP], output_data["Region-Loss"])
        g = seaborn.jointplot(data=output_data, x=MSP, y="Region-Loss", hue="Stage", hue_order=stage_list, palette=step00.precancer_color_code, height=18, ratio=5, kind="scatter")
        if p < 1e-3:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        g.set_axis_labels(MSP, "Number of somatic CNV-Loss region (Gbp)")
        figures.append(f"Joint_Precancer-Loss_{MSP}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

        r, p = scipy.stats.spearmanr(output_data[MSP], output_data["Region-Gain"])
        g = seaborn.jointplot(data=output_data, x=MSP, y="Region-Gain", hue="Stage", hue_order=stage_list, palette=step00.precancer_color_code, height=18, ratio=5, kind="scatter")
        if p < 1e-3:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        g.set_axis_labels(MSP, "Number of somatic CNV-Gain region (Gbp)")
        figures.append(f"Joint_Precancer-Gain_{MSP}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

    output_data = output_data.loc[(output_data["Stage"] == "Precancer")]
    for MSP in tqdm.tqdm(step00.sharing_columns[:2]):
        output_data[MSP] = list(map(lambda x: clinical_data.loc[x, MSP], output_data["Patient"]))

        r, p = scipy.stats.spearmanr(output_data[MSP], output_data["Region"])
        g = seaborn.jointplot(data=output_data, x=MSP, y="Region", color=step00.precancer_color_code["Precancer"], height=18, ratio=5, kind="reg")
        if p < 1e-3:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        g.set_axis_labels(MSP, "Number of somatic CNV region (Gbp)")
        figures.append(f"Joint_PrecancerOnly_{MSP}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.regplot(data=output_data, x=MSP, y="Region", fit_reg=True, scatter=True, color=step00.precancer_color_code["Precancer"], truncate=False, ax=ax)

        if p < 1e-3:
            matplotlib.pyplot.text(get_middle(output_data[MSP]), get_middle(output_data["Region"]), f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            matplotlib.pyplot.text(get_middle(output_data[MSP]), get_middle(output_data["Region"]), f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.xlabel("PSM")
        matplotlib.pyplot.ylabel("Number of somatic CNV region (Gbp)")
        matplotlib.pyplot.tight_layout()

        figures.append(f"Scatter_PrecancerOnly_{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

        r, p = scipy.stats.spearmanr(output_data[MSP], output_data["Region-Loss"])
        g = seaborn.jointplot(data=output_data, x=MSP, y="Region-Loss", color=step00.precancer_color_code["Precancer"], height=18, ratio=5, kind="reg")
        if p < 1e-3:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        g.set_axis_labels(MSP, "Number of somatic CNV-Loss region (Gbp)")
        figures.append(f"Joint_PrecancerOnly-Loss_{MSP}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.regplot(data=output_data, x=MSP, y="Region-Loss", fit_reg=True, scatter=True, color=step00.precancer_color_code["Precancer"], truncate=False, ax=ax)

        if p < 1e-3:
            matplotlib.pyplot.text(get_middle(output_data[MSP]), get_middle(output_data["Region-Loss"]), f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            matplotlib.pyplot.text(get_middle(output_data[MSP]), get_middle(output_data["Region-Loss"]), f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.xlabel("PSM")
        matplotlib.pyplot.ylabel("Number of somatic CNV-Loss region (Gbp)")
        matplotlib.pyplot.tight_layout()

        figures.append(f"Scatter_PrecancerOnly-Loss_{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

        r, p = scipy.stats.spearmanr(output_data[MSP], output_data["Region-Gain"])
        g = seaborn.jointplot(data=output_data, x=MSP, y="Region-Gain", color=step00.precancer_color_code["Precancer"], height=18, ratio=5, kind="reg")

        if p < 1e-3:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            g.fig.text(0.5, 0.5, f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        g.set_axis_labels(MSP, "Number of somatic CNV-Gain region (Gbp)")
        figures.append(f"Joint_PrecancerOnly-Gain_{MSP}.pdf")
        g.savefig(figures[-1])
        matplotlib.pyplot.close(g.fig)

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.regplot(data=output_data, x=MSP, y="Region-Gain", fit_reg=True, scatter=True, color=step00.precancer_color_code["Precancer"], truncate=False, ax=ax)

        if p < 1e-3:
            matplotlib.pyplot.text(get_middle(output_data[MSP]), get_middle(output_data["Region-Gain"]), f"r={r:.3f}, p={p:.1e}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})
        else:
            matplotlib.pyplot.text(get_middle(output_data[MSP]), get_middle(output_data["Region-Gain"]), f"r={r:.3f}, p={p:.3f}", color="k", fontsize="small", horizontalalignment="center", verticalalignment="center", bbox={"alpha": 0.3, "color": "white"})

        matplotlib.pyplot.xlabel("PSM")
        matplotlib.pyplot.ylabel("Number of somatic CNV-Gain region (Gbp)")
        matplotlib.pyplot.tight_layout()

        figures.append(f"Scatter_PrecancerOnly-Gain_{MSP}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
