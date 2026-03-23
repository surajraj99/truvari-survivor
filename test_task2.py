import sys
import os

# Add truvari to path
sys.path.insert(0, os.path.abspath("."))

from truvari.survivor import parse_args, get_callers, tagged_variant_stream, Caller
import truvari

def test_parsing():
    args = ["-i", "input1.vcf.gz", "input2.vcf.gz", "-o", "output.vcf", "-r", "1000"]
    parsed = parse_args(args)
    assert parsed.input == ["input1.vcf.gz", "input2.vcf.gz"]
    assert parsed.output == "output.vcf"
    assert parsed.refdist == 1000
    print("test_parsing passed")

def test_get_callers():
    input_vcfs = ["input1.vcf.gz", "input2.vcf.gz"]
    callers = get_callers(input_vcfs, None)
    assert len(callers) == 2
    assert callers[0].name == "input1.vcf.gz"
    assert callers[0].index == 0
    assert callers[1].index == 1
    print("test_get_callers passed")

def test_list_parsing():
    with open("vcf_list.txt", "w") as f:
        f.write("input3.vcf.gz\n")
    
    callers = get_callers(["input1.vcf.gz"], "vcf_list.txt")
    assert len(callers) == 2
    assert callers[1].name == "input3.vcf.gz"
    assert callers[1].index == 1
    os.remove("vcf_list.txt")
    print("test_list_parsing passed")

def test_stream():
    vcf1 = "repo_utils/test_files/variants/input1.vcf.gz"
    vcf2 = "repo_utils/test_files/variants/input2.vcf.gz"
    callers = [
        Caller(name="vcf1", vcf_path=vcf1, index=0),
        Caller(name="vcf2", vcf_path=vcf2, index=1)
    ]
    params = truvari.VariantParams()
    count = 0
    for record in tagged_variant_stream(callers, params):
        assert hasattr(record, "caller_id")
        assert record.caller_id in [0, 1]
        count += 1
    assert count > 0
    print(f"test_stream passed with {count} records")

if __name__ == "__main__":
    try:
        test_parsing()
        test_get_callers()
        test_list_parsing()
        test_stream()
        print("All tests passed!")
    except Exception as e:
        print(f"Test failed: {e}")
        sys.exit(1)
