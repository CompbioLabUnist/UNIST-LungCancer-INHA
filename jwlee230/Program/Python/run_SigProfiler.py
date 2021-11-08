"""
run_SigProfiler.py: run SigProfilerExtractor as https://github.com/AlexandrovLab/SigProfilerExtractor
"""
import argparse
import os
import os.path
import tempfile
from SigProfilerExtractor import sigpro as sig
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input input Mutect2 VCF file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output directory", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".vcf"), args.input)):
        raise ValueError("Input files must end with .vcf!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be a positive integer!!")

    with tempfile.TemporaryDirectory(dir=step00.tmpfs) as tmpdir:
        input_directory = os.path.join(tmpdir, "input")
        os.makedirs(input_directory, mode=775)

        for input_file in args.input:
            os.system("cp -v {0} {1}".format(input_file, os.path.join(input_directory, step00.get_id(input_file) + ".vcf")))

        sig.sigProfilerExtractor("vcf", args.output, input_data=tmpdir, reference_genome="GRCh38", opportunity_genome="GRCh38", exome=True, cpu=args.cpus)
