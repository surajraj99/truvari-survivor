"""
N-way VCF merging with SURVIVOR-style support tracking
"""
import argparse
import truvari

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="survivor", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", nargs="+", help="Input VCFs")
    parser.add_argument("-l", "--list", type=str, help="File list of VCFs")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output VCF")
    return parser.parse_args(args)

def survivor_main(args):
    """
    Main entrypoint for survivor
    """
    args = parse_args(args)
    print("Truvari Survivor Placeholder")
