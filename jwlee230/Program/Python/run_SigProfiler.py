"""
run_SigProfiler.py: run SigProfilerExtractor as https://github.com/AlexandrovLab/SigProfilerExtractor
"""
import argparse
from SigProfilerExtractor import sigpro as sig

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input input Mutect2 VCF dirctory", type=str)
    parser.add_argument("output", help="Output directory", type=str)
    parser.add_argument("--cpus", help="CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if args.cpus < 1:
        raise ValueError("CPUs must be a positive integer!!")

    sig.sigProfilerExtractor("vcf", args.output, input_data=args.input, reference_genome="GRCh38", opportunity_genome="GRCh38", exome=True, cpu=args.cpus)
