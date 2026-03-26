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
    thresg.add_argument("-p", "--pctseq", type=truvari.restricted_float, default=0.0,
                        help="Min percent sequence similarity. Set to 0 to ignore. (%(default)s)")
    thresg.add_argument("-P", "--pctsize", type=truvari.restricted_float, default=0.50,
                        help="Min pct allele size similarity (minvarsize/maxvarsize) (%(default)s)")
    thresg.add_argument("-O", "--pctovl", type=truvari.restricted_float, default=0.50,
                        help="Min pct reciprocal overlap (%(default)s)")
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
    thresg.add_argument("-B", "--bnddist", type=int, default=None,
                        help="Maximum distance allowed between BNDs (refdist; -1=off)")
    thresg.add_argument("--strandignore", action="store_true", default=True,
                        help="Ignore strand mismatches during comparison (%(default)s)")
    thresg.add_argument("--no-strandignore", action="store_false", dest="strandignore",
                        help="Do not ignore strand mismatches during comparison")
    thresg.add_argument("--dynthresh", type=str, default=None,
                        help="Dynamic thres params (min_diff,max_diff,s_min,s_max) overrides pctsize/pctseq (off)")

    filteg = parser.add_argument_group("Filtering Arguments")
    filteg.add_argument("-s", "--sizemin", type=truvari.restricted_int, default=50,
                        help="Minimum variant size to consider for comparison (%(default)s)")
    filteg.add_argument("-S", "--sizemax", type=int, default=50000,
                        help="Maximum variant size to consider for comparison (-1 = off; %(default)s)")
    parser.add_argument("--passonly", action="store_true", default=False,
                        help="Only consider calls with FILTER == PASS")

    parser.add_argument("--debug", action="store_true", default=False,
                        help="Verbose logging")
    parser.add_argument("--log", type=str, default=None,
                        help="Log file to save debug output")

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
    header.add_line('##INFO=<ID=CALLER_IDS,Number=.,Type=String,Description="Original variant IDs from each supporting caller">')
    
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
    if args.bnddist is None:
        args.bnddist = args.refdist
    
    # Setup console logging. If --debug is True, this sets root to DEBUG.
    # Otherwise, it sets root to INFO.
    truvari.setup_logging(args.debug, show_version=True)
    
    # If we have a log file, we need to ensure the root logger is at DEBUG 
    # so that debug messages are generated at all.
    if args.log:
        if not args.debug:
            # Bump root to DEBUG to catch everything for the file
            logging.getLogger().setLevel(logging.DEBUG)
            # But explicitly restrict console handlers to INFO
            for handler in logging.getLogger().handlers:
                if isinstance(handler, logging.StreamHandler):
                    handler.setLevel(logging.INFO)
        
        # Add the FileHandler at DEBUG level
        fh = logging.FileHandler(args.log)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
        logging.getLogger().addHandler(fh)

    if not args.input and not args.list:
        logging.error("Must specify either --input or --list")
        sys.exit(100)

    callers = get_callers(args.input, args.list)

    params = truvari.VariantParams(args=args,
                                   short_circuit=False,
                                   skip_gt=True,
                                   sizefilt=args.sizemin)
    reference = pysam.FastaFile(params.reference) if params.reference else None

    # Extra attributes needed
    params.sorter = SORTS['first']
    params.gt = 'off'
    params.chain = True # Enable chaining by default for survivor-style merging
    params.chunksize = 10000000 # Very large chunksize to prevent splitting variants into different chunks

    out_header = create_survivor_header(callers)
    out_vcf = pysam.VariantFile(args.output, 'w', header=out_header)

    # Use file_zipper to get variants from all callers in a stream
    file_iters = [('base', tagged_variant_stream_single(caller, params)) for caller in callers]
    variant_stream = truvari.file_zipper(*file_iters)

    # Collect all variants from all callers and group by chromosome
    # This prevents split-chromosome issues due to non-lexicographical sorting in VCFs
    all_chrom_variants = defaultdict(list)
    for key, var in variant_stream:
        all_chrom_variants[var.chrom].append(var)

    for chrom, variants in all_chrom_variants.items():
        logging.info("Processing chromosome %s", chrom)
        
        # Ensure they are sorted by position
        variants.sort(key=lambda x: x.pos)
        
        used = [False] * len(variants)
        group_count = 0
        
        for i, var in enumerate(variants):
            if used[i]:
                continue
            
            # Start a new group
            group_count += 1
            current_group = [var]
            used[i] = True
            
            # Find all matching variants (chaining)
            changed = True
            while changed:
                changed = False
                for j, candidate in enumerate(variants):
                    if used[j]:
                        continue
                    
                    # Check if candidate matches ANY member of current_group
                    matched = False
                    for member in current_group:
                        vstart, vend = member.boundaries()
                        cstart, cend = candidate.boundaries()
                        
                        # Within range check
                        if truvari.overlaps(vstart - params.refdist, vend + params.refdist, cstart, cend):
                            mat = member.match(candidate)
                            if mat.state:
                                matched = True
                                break
                    
                    if matched:
                        current_group.append(candidate)
                        used[j] = True
                        changed = True
            
            if args.debug:
                logging.debug("Group %d has %d members: %s", group_count, len(current_group), 
                              ", ".join([f"{v!r}" for v in current_group]))
            
            # current_group now has all linked variants
            # Calculate median POS and END
            import statistics
            import re
            
            # Median position and stop
            med_pos = int(statistics.median([v.pos for v in current_group]))
            med_stop = int(statistics.median([v.stop for v in current_group]))
            
            # Determine the consensus SVTYPE (most frequent normalized type)
            from collections import Counter
            type_counts = Counter([v.var_type().name for v in current_group])
            consensus_type = type_counts.most_common(1)[0][0]

            # For BND/TRA, also find consensus CHR2 and second positions
            consensus_chr2 = None
            mate_positions = []
            if consensus_type in ["BND", "TRA"]:
                chr2_counts = Counter()
                for v in current_group:
                    try:
                        mc, mp = v.bnd_position()
                        chr2_counts[mc] += 1
                        mate_positions.append(mp)
                    except ValueError:
                        continue
                if chr2_counts:
                    consensus_chr2 = chr2_counts.most_common(1)[0][0]

            # We pick the representative as the one closest to the median
            rep = min(current_group, key=lambda x: abs(x.pos - med_pos))
            
            # Create a new record in the output
            new_record = out_vcf.new_record()
            new_record.chrom = rep.chrom
            new_record.pos = med_pos
            new_record.id = rep.id
            
            # Ensure REF consistency
            if reference:
                try:
                    new_record.ref = reference.fetch(rep.chrom, med_pos, med_pos + 1)
                except (ValueError, KeyError, IndexError):
                    new_record.ref = rep.ref
            else:
                new_record.ref = rep.ref
            
            # Handle BND ALTs
            if consensus_type == "BND":
                if mate_positions:
                    med_mate_pos = int(statistics.median(mate_positions))
                    # Reconstruct ALT using med_mate_pos and consensus_chr2
                    # We take the bracket style from the representative
                    new_alt = re.sub(r'([\[\]])[^\[\]]+:[0-9]+([\[\]])', 
                                     fr'\g<1>{consensus_chr2}:{med_mate_pos}\g<2>', 
                                     rep.alts[0])
                    new_record.alts = (new_alt,)
                else:
                    new_record.alts = rep.alts
            else:
                new_record.alts = rep.alts

            new_record.qual = rep.qual
            
            # Copy filters
            for f in rep.filter.keys():
                new_record.filter.add(f)
            
            # Copy INFO from representative
            for k, v in rep.info.items():
                if k in out_header.info and k not in ["END", "SVLEN", "SVTYPE", "CHR2"]:
                    try:
                        new_record.info[k] = v
                    except (TypeError, ValueError):
                        continue
            
            # Set consensus/calculated fields
            new_record.info["SVTYPE"] = consensus_type
            if consensus_type == "TRA":
                new_record.info["CHR2"] = consensus_chr2
                if mate_positions:
                    new_record.stop = int(statistics.median(mate_positions))
                else:
                    new_record.stop = med_stop
            else:
                new_record.stop = med_stop
            
            # Calculate SVLEN
            # For standard SVs, SVLEN is often end - pos (negative for DEL)
            # Truvari's var_size() is always positive. 
            # We'll follow the convention of the consensus type if possible
            if consensus_type in ["BND", "TRA"]:
                svlen = 0
            else:
                svlen = med_stop - med_pos
                if consensus_type == "DEL":
                    svlen = -svlen
                elif consensus_type == "INS":
                    # For insertions, end is usually pos + 1, so med_stop - med_pos = 1
                    # We should take the median of the actual SVLENs if they exist
                    sizes = []
                    for v in current_group:
                        sizes.append(v.var_size())
                    svlen = int(statistics.median(sizes))
            
            new_record.info["SVLEN"] = svlen

            # Calculate SURVIVOR fields
            num_callers = len(callers)
            supp_vec = [0] * num_callers
            caller_ids = set()
            caller_variant_ids = []
            seen_caller_ids = set()
            
            for v in current_group:
                if hasattr(v, 'caller_id'):
                    caller_ids.add(v.caller_id)
                    supp_vec[v.caller_id] = 1
                    # Store original variant ID if it exists
                    if v.id and v.id != '.':
                        caller_variant_ids.append(f"{callers[v.caller_id].name}:{v.id}")
            
            new_record.info["SUPP"] = len(caller_ids)
            new_record.info["SUPP_VEC"] = "".join(map(str, supp_vec))
            new_record.info["CALLERS"] = [callers[idx].name for idx in sorted(list(caller_ids))]
            if caller_variant_ids:
                new_record.info["CALLER_IDS"] = caller_variant_ids
            
            # Set genotypes for each caller
            caller_to_record = {}
            for v in current_group:
                if hasattr(v, 'caller_id'):
                    if v.caller_id not in caller_to_record:
                        caller_to_record[v.caller_id] = v
            
            for caller in callers:
                if caller.index in caller_to_record:
                    src = caller_to_record[caller.index]
                    try:
                        new_record.samples[caller.name]['GT'] = src.samples[0]['GT']
                    except (KeyError, IndexError):
                        new_record.samples[caller.name]['GT'] = (None, None)
                else:
                    new_record.samples[caller.name]['GT'] = (None, None)
            
            out_vcf.write(new_record)

    out_vcf.close()
    logging.info("Finished survivor")
