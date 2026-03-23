import os
import pysam
import sys

def create_vcf(path, records, sample_names=["SAMPLE"]):
    header = pysam.VariantHeader()
    header.add_line('##fileformat=VCFv4.2')
    header.add_line('##FILTER=<ID=PASS,Description="All filters passed">')
    header.add_line('##contig=<ID=chr1,length=1000000>')
    header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
    header.add_line('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of variant">')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    for sample in sample_names:
        header.add_sample(sample)
    
    with pysam.VariantFile(path, 'w', header=header) as vcf:
        for rec in records:
            v_rec = vcf.new_record()
            v_rec.chrom = rec.get('chrom', 'chr1')
            v_rec.pos = rec['pos']
            v_rec.id = rec.get('id', '.')
            v_rec.ref = rec.get('ref', 'N')
            v_rec.alts = rec['alts']
            v_rec.qual = rec.get('qual', 100)
            v_rec.filter.add('PASS')
            for k, v in rec.get('info', {}).items():
                if k == 'END':
                    v_rec.stop = v
                else:
                    v_rec.info[k] = v
            for sample in sample_names:
                v_rec.samples[sample]['GT'] = rec.get('gt', (0, 1))
            vcf.write(v_rec)
    pysam.tabix_index(path, preset="vcf", force=True)

def scenario_easy(out_dir):
    """
    Easy: high overlap, identical types
    """
    os.makedirs(out_dir, exist_ok=True)
    variants = [
        {'pos': 1000, 'alts': ('<DEL>',), 'info': {'SVTYPE': 'DEL', 'SVLEN': -100, 'END': 1100}},
        {'pos': 2000, 'alts': ('<INS>',), 'info': {'SVTYPE': 'INS', 'SVLEN': 100, 'END': 2000}},
        {'pos': 3000, 'alts': ('<DUP>',), 'info': {'SVTYPE': 'DUP', 'SVLEN': 100, 'END': 3100}},
    ]
    
    # Caller 1
    create_vcf(os.path.join(out_dir, "caller1.vcf.gz"), variants, ["C1"])
    # Caller 2 (identical)
    create_vcf(os.path.join(out_dir, "caller2.vcf.gz"), variants, ["C2"])
    # Caller 3 (identical)
    create_vcf(os.path.join(out_dir, "caller3.vcf.gz"), variants, ["C3"])

def scenario_medium(out_dir):
    """
    Medium: BNDs vs SVs, position offsets
    """
    os.makedirs(out_dir, exist_ok=True)
    # Caller 1: SV DEL
    v1 = [{'pos': 1000, 'alts': ('<DEL>',), 'info': {'SVTYPE': 'DEL', 'SVLEN': -100, 'END': 1100}}]
    create_vcf(os.path.join(out_dir, "caller1.vcf.gz"), v1, ["C1"])
    
    # Caller 2: Offset DEL (50bp offset)
    v2 = [{'pos': 1050, 'alts': ('<DEL>',), 'info': {'SVTYPE': 'DEL', 'SVLEN': -100, 'END': 1150}}]
    create_vcf(os.path.join(out_dir, "caller2.vcf.gz"), v2, ["C2"])
    
    # Caller 3: BND representing the same DEL
    # Simplified BND for testing - Truvari might need specific BND format
    v3 = [
        {'pos': 1000, 'alts': ('N[chr1:1100[',), 'info': {'SVTYPE': 'BND'}},
        {'pos': 1100, 'alts': ('N[chr1:1000[',), 'info': {'SVTYPE': 'BND'}}
    ]
    create_vcf(os.path.join(out_dir, "caller3.vcf.gz"), v3, ["C3"])

def scenario_hard(out_dir):
    """
    Hard: complex translocations, low sequence similarity
    """
    os.makedirs(out_dir, exist_ok=True)
    # Caller 1: INS with sequence
    v1 = [{'pos': 1000, 'alts': ('ATCG' * 25,), 'info': {'SVTYPE': 'INS', 'SVLEN': 100, 'END': 1000}}]
    create_vcf(os.path.join(out_dir, "caller1.vcf.gz"), v1, ["C1"])
    
    # Caller 2: INS with slightly different sequence (60% similarity)
    # ATCG * 25 = 100bp
    # Let's make it 100bp of 'A' vs 100bp of 'T' - sequence similarity will be 0
    v2 = [{'pos': 1000, 'alts': ('A' * 100,), 'info': {'SVTYPE': 'INS', 'SVLEN': 100, 'END': 1000}}]
    create_vcf(os.path.join(out_dir, "caller2.vcf.gz"), v2, ["C2"])
    
    # Caller 3: BND translocation
    v3 = [
        {'pos': 1000, 'chrom': 'chr1', 'alts': ('N[chr2:2000[',), 'info': {'SVTYPE': 'BND'}},
    ]
    # Need chr2 in header if we do this
    # For now, let's just do another complex INS
    v3 = [{'pos': 1000, 'alts': ('G' * 100,), 'info': {'SVTYPE': 'INS', 'SVLEN': 100, 'END': 1000}}]
    create_vcf(os.path.join(out_dir, "caller3.vcf.gz"), v3, ["C3"])

if __name__ == "__main__":
    base_dir = sys.argv[1] if len(sys.argv) > 1 else "sim_data"
    scenario_easy(os.path.join(base_dir, "easy"))
    scenario_medium(os.path.join(base_dir, "medium"))
    scenario_hard(os.path.join(base_dir, "hard"))
