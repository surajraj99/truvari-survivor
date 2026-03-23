"""
N-way VCF merging with SURVIVOR-style support tracking
"""
import os
import sys
import logging
import argparse
from dataclasses import dataclass
from collections import defaultdict

import pysam
import truvari
from truvari.collapse import collapse_chunk, SORTS

@dataclass
class Caller:
    """
    Source VCF information
    """
    name: str
    vcf_path: str
    index: int

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="survivor", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", nargs="+", help="Input VCFs")
    parser.add_argument("-l", "--list", type=str, help="File list of VCFs")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output VCF")

    thresg = parser.add_argument_group("Comparison Threshold Arguments")
    thresg.add_argument("-r", "--refdist", type=truvari.restricted_int, default=500,
                        help="Max reference location distance (%(default)s)")
    thresg.add_argument("-p", "--pctseq", type=truvari.restricted_float, default=0.70,
                        help="Min percent sequence similarity. Set to 0 to ignore. (%(default)s)")
    thresg.add_argument("-P", "--pctsize", type=truvari.restricted_float, default=0.70,
                        help="Min pct allele size similarity (minvarsize/maxvarsize) (%(default)s)")
    thresg.add_argument("-O", "--pctovl", type=truvari.restricted_float, default=0.0,
                        help="Min pct reciprocal overlap (%(default)s) for DEL events")
    thresg.add_argument("-t", "--typeignore", action="store_true", default=False,
                        help="Variant types don't need to match to compare (%(default)s)")
    thresg.add_argument("-n", "--no-roll", action="store_true",
                        help="Turn off rolling sequence similarity")
    thresg.add_argument("-m", "--max-resolve", type=truvari.restricted_int, default=25000,
                        help="Maximum size of variant to attempt to sequence resolve ($(default)s)")
    thresg.add_argument("-D", "--decompose", action="store_true",
                        help="Allow decomposition for SV to BND comparison (%(default)s)")
    thresg.add_argument("-d", "--dup-to-ins", action="store_true",
                        help="Assume DUP svtypes are INS (%(default)s)")
    thresg.add_argument("-B", "--bnddist", type=int, default=100,
                        help="Maximum distance allowed between BNDs (%(default)s; -1=off)")
    thresg.add_argument("--dynthresh", type=str, default=None,
                        help="Dynamic thres params (min_diff,max_diff,s_min,s_max) overrides pctsize/pctseq (off)")

    filteg = parser.add_argument_group("Filtering Arguments")
    filteg.add_argument("-s", "--sizemin", type=truvari.restricted_int, default=50,
                        help="Minimum variant size to consider for comparison (%(default)s)")
    filteg.add_argument("-S", "--sizemax", type=int, default=50000,
                        help="Maximum variant size to consider for comparison (-1 = off; %(default)s)")
    filteg.add_argument("--passonly", action="store_true", default=False,
                        help="Only consider calls with FILTER == PASS")

    parser.add_argument("--debug", action="store_true", default=False,
                        help="Verbose logging")

    args = parser.parse_args(args)
    if args.dynthresh is not None:
        args.dynthresh = list(map(int, args.dynthresh.split(',')))
        tmp = truvari.VariantParams(dynthresh=args.dynthresh).calc_dyn_thresh(args.sizemin)
        args.pctseq = tmp
        args.pctsize = tmp
    return args

def get_callers(input_vcfs, vcf_list):
    """
    Load callers from input and/or list
    """
    callers = []
    if input_vcfs:
        for i, vcf in enumerate(input_vcfs):
            callers.append(Caller(name=vcf, vcf_path=vcf, index=i))

    if vcf_list:
        with open(vcf_list, 'r') as fh:
            for line in fh:
                vcf = line.strip()
                if vcf:
                    callers.append(Caller(name=vcf, vcf_path=vcf, index=len(callers)))
    return callers

def tagged_variant_stream(callers, params):
    """
    Iterate over callers and yield tagged records
    """
    for caller in callers:
        with truvari.VariantFile(caller.vcf_path, params=params) as vcf:
            for record in vcf:
                record.caller_id = caller.index
                yield record

def tagged_variant_stream_single(caller, params):
    """
    Yields records from a single caller, tagged with caller_id
    """
    with truvari.VariantFile(caller.vcf_path, params=params) as vcf:
        for record in vcf:
            record.caller_id = caller.index
            yield record

def create_survivor_header(callers):
    """
    Create a VCF header for SURVIVOR output by merging all input headers
    """
    header = pysam.VariantHeader()
    
    # Add standard SURVIVOR fields first so they have our descriptions
    header.add_line('##INFO=<ID=SUPP,Number=1,Type=Integer,Description="Number of callers supporting the variant">')
    header.add_line('##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Vector of callers supporting the variant">')
    header.add_line('##INFO=<ID=CALLERS,Number=.,Type=String,Description="Names of callers supporting the variant">')
    
    # Add standard GT field
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    
    # Add standard SV fields that are often expected
    header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
    header.add_line('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">')

    for caller in callers:
        with truvari.VariantFile(caller.vcf_path) as vcf:
            for record in vcf.header.records:
                if record.type in ["CONTIG", "FILTER", "FORMAT", "INFO"]:
                    # Skip if already exists to avoid duplicates or overwriting our definitions
                    rid = record.get("ID")
                    if rid is None:
                        continue
                    if record.type == "INFO" and rid in header.info:
                        continue
                    if record.type == "FORMAT" and rid in header.formats:
                        continue
                    if record.type == "FILTER" and rid in header.filters:
                        continue
                    if record.type == "CONTIG" and rid in header.contigs:
                        continue
                    header.add_record(record)
        header.add_sample(caller.name)
        
    return header

def annotate_survivor(new_record, group, callers):
    """
    Add SURVIVOR fields to the entry
    """
    num_callers = len(callers)
    supp_vec = [0] * num_callers
    
    caller_ids = set()
    if hasattr(group.entry, 'caller_id'):
        caller_ids.add(group.entry.caller_id)
        supp_vec[group.entry.caller_id] = 1
    
    for match in group.matches:
        if hasattr(match.comp, 'caller_id'):
            caller_ids.add(match.comp.caller_id)
            supp_vec[match.comp.caller_id] = 1
            
    new_record.info["SUPP"] = len(caller_ids)
    new_record.info["SUPP_VEC"] = "".join(map(str, supp_vec))
    new_record.info["CALLERS"] = [callers[i].name for i in sorted(list(caller_ids))]

def survivor_main(args):
    """
    Main entrypoint for survivor
    """
    args = parse_args(args)
    truvari.setup_logging(args.debug, show_version=True)

    if not args.input and not args.list:
        logging.error("Must specify either --input or --list")
        sys.exit(100)

    callers = get_callers(args.input, args.list)

    params = truvari.VariantParams(args=args,
                                   short_circuit=True,
                                   skip_gt=True,
                                   sizefilt=args.sizemin)
    # Extra attributes needed for collapse_chunk
    params.sorter = SORTS['first']
    params.gt = 'off'
    params.chain = False
    params.hap = False
    params.consolidate = False

    out_header = create_survivor_header(callers)
    out_vcf = pysam.VariantFile(args.output, 'w', header=out_header)

    # Use chunker to get variants from all callers
    file_iters = [('base', tagged_variant_stream_single(caller, params)) for caller in callers]
    chunks = truvari.chunker(params, *file_iters)

    for chunk, chunk_id in chunks:
        # Since we use 'base' for all, they are all in chunk['base']
        collapsed_groups = collapse_chunk((chunk, chunk_id), params)
        
        for group in collapsed_groups:
            # group is a CollapsedCalls object
            if group.entry is None:
                continue
                
            # Create a new record in the output
            new_record = out_vcf.new_record()
            new_record.chrom = group.entry.chrom
            new_record.pos = group.entry.pos
            new_record.id = group.entry.id
            new_record.ref = group.entry.ref
            new_record.alts = group.entry.alts
            new_record.qual = group.entry.qual
            
            # Copy filters
            for f in group.entry.filter.keys():
                new_record.filter.add(f)
            
            # Copy INFO
            for k, v in group.entry.info.items():
                if k in out_header.info and k != "END":
                    try:
                        new_record.info[k] = v
                    except (TypeError, ValueError):
                        continue
            
            # Properly set END/stop
            new_record.stop = group.entry.stop

            # Annotate with SURVIVOR fields
            annotate_survivor(new_record, group, callers)
            
            # Set genotypes for each caller
            caller_to_record = {}
            if hasattr(group.entry, 'caller_id'):
                caller_to_record[group.entry.caller_id] = group.entry
            for match in group.matches:
                if hasattr(match.comp, 'caller_id'):
                    if match.comp.caller_id not in caller_to_record:
                        caller_to_record[match.comp.caller_id] = match.comp
            
            for caller in callers:
                if caller.index in caller_to_record:
                    src = caller_to_record[caller.index]
                    # SURVIVOR format: just the first sample's GT
                    try:
                        new_record.samples[caller.name]['GT'] = src.samples[0]['GT']
                    except (KeyError, IndexError):
                        new_record.samples[caller.name]['GT'] = (None, None)
                else:
                    new_record.samples[caller.name]['GT'] = (None, None)
            
            out_vcf.write(new_record)

    out_vcf.close()
    logging.info("Finished survivor")
