"""
calculate_CNV_MSP.py: calculate correlation between CNV DEG and MSP
"""
import argparse
import itertools
import multiprocessing
import typing
import numpy
import pandas
import scipy.stats
import tqdm
import step00

input_data = pandas.DataFrame()


def correlation_stage(stage: str, MSP: str, gene: str) -> typing.Tuple[float, float, float, float, float, float]:
    tmp_data = input_data.loc[(input_data["Stage"] == stage), [MSP, gene]]
    if tmp_data.shape[0] < 3:
        return 0.0, 0.0, 0.0, 1.0, 0.0, 0.0
    elif (numpy.std(tmp_data[MSP]) == 0.0) or (numpy.std(tmp_data[gene]) == 0.0):
        return 0.0, 0.0, 0.0, 1.0, 0.0, 0.0
    return scipy.stats.linregress(tmp_data[MSP], tmp_data[gene])


def correlation_precancer(MSP: str, gene: str) -> typing.Tuple[float, float, float, float, float, float]:
    tmp_data = input_data.loc[~(input_data["Stage"].isin({"Normal", "Primary"})), [MSP, gene]]
    if tmp_data.shape[0] < 3:
        return 0.0, 0.0, 0.0, 1.0, 0.0, 0.0
    elif (numpy.std(tmp_data[MSP]) == 0.0) or (numpy.std(tmp_data[gene]) == 0.0):
        return 0.0, 0.0, 0.0, 1.0, 0.0, 0.0
    return scipy.stats.linregress(tmp_data[MSP], tmp_data[gene])


def correlation_all(MSP: str, gene: str) -> typing.Tuple[float, float, float, float, float, float]:
    tmp_data = input_data.loc[:, [MSP, gene]]
    if tmp_data.shape[0] < 3:
        return 0.0, 0.0, 0.0, 1.0, 0.0, 0.0
    elif (numpy.std(tmp_data[MSP]) == 0.0) or (numpy.std(tmp_data[gene]) == 0.0):
        return 0.0, 0.0, 0.0, 1.0, 0.0, 0.0
    return scipy.stats.linregress(tmp_data[MSP], tmp_data[gene])
