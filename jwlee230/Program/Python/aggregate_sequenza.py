"""
aggregate_sequenza.py: aggregate sequenza results
"""
import argparse
import os
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input sequenza TAR.gz files", type=str, nargs="+")
    parser.add_argument("output", help="Output PNG file", type=str)
    parser.add_argument("--tmpfs", help="Path to temporary directory", type=str, default="/tmpfs")

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tar.gz"), args.input)):
        raise ValueError("Input must end with .tar.gz!!")
    elif not args.output.endswith(".png"):
        raise ValueError("Output must end with .png!!")

    args.input.sort()
    print(len(args.input))
