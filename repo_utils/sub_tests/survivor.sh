# ------------------------------------------------------------
#                                 survivor
# ------------------------------------------------------------
run survivor_help $truv survivor -h
assert_in_stdout "N-way VCF merging"

# Test input parsing with both -i and -l
VCF1=$INDIR/variants/input1.vcf.gz
VCF2=$INDIR/variants/input2.vcf.gz
VCF3=$INDIR/variants/input3.vcf.gz
echo $VCF3 > vcf_list.txt

run survivor_input $truv survivor -i $VCF1 $VCF2 -l vcf_list.txt -o $OD/out.vcf
assert_in_stderr "Finished survivor"
assert_exit_code 0

rm vcf_list.txt

# Setup simulation data
run survivor_sim python3 repo_utils/test_files/scripts/sim_survivor.py $OD/sim_data
assert_exit_code 0

# 1. Easy Scenario: Identical variants across 3 callers
VCFS="$OD/sim_data/easy/caller1.vcf.gz $OD/sim_data/easy/caller2.vcf.gz $OD/sim_data/easy/caller3.vcf.gz"
run survivor_easy $truv survivor -i $VCFS -o $OD/easy.vcf
assert_exit_code 0
# Check SUPP=3 for all records and verify SUPP_VEC and CALLERS
run survivor_easy_check python3 -c "
import pysam
v = pysam.VariantFile('$OD/easy.vcf')
all_good = True
for r in v:
    if r.info['SUPP'] != 3:
        all_good = False
    if r.info['SUPP_VEC'] != '111':
        all_good = False
    if len(r.info['CALLERS']) != 3:
        all_good = False
print(all_good)
"
assert_stdout "True"

# 2. Medium Scenario: SV to BND with position offsets
VCFS="$OD/sim_data/medium/caller1.vcf.gz $OD/sim_data/medium/caller2.vcf.gz $OD/sim_data/medium/caller3.vcf.gz"
# Need -D for SV-to-BND comparison
run survivor_medium $truv survivor -i $VCFS -o $OD/medium.vcf -D
assert_exit_code 0
# They should cluster into one or more groups.
# C1 and C2 DELs should match. C3 BNDs should match DELs with -D.
run survivor_medium_check python3 -c "
import pysam
v = pysam.VariantFile('$OD/medium.vcf')
found_supp3 = False
for r in v:
    if r.info['SUPP'] == 3:
        found_supp3 = True
print(found_supp3)
"
assert_stdout "True"

# 3. Hard Scenario: Different sequences for INS
VCFS="$OD/sim_data/hard/caller1.vcf.gz $OD/sim_data/hard/caller2.vcf.gz $OD/sim_data/hard/caller3.vcf.gz"
# With default pctseq=0.70, these should NOT match (seq similarity is 0)
run survivor_hard $truv survivor -i $VCFS -o $OD/hard.vcf
assert_exit_code 0
# Should have 3 records with SUPP=1
run survivor_hard_check python3 -c "
import pysam
v = pysam.VariantFile('$OD/hard.vcf')
count = sum(1 for r in v if r.info['SUPP'] == 1)
print(count)
"
assert_stdout "3"

# Try hard again with pctseq=0 to allow merging regardless of sequence
run survivor_hard_no_seq $truv survivor -i $VCFS -o $OD/hard_no_seq.vcf -p 0
assert_exit_code 0
# Should now have 1 record with SUPP=3
run survivor_hard_no_seq_check python3 -c "
import pysam
v = pysam.VariantFile('$OD/hard_no_seq.vcf')
found_supp3 = False
for r in v:
    if r.info['SUPP'] == 3:
        found_supp3 = True
print(found_supp3)
"
assert_stdout "True"
