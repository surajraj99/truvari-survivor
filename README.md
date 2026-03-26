[![PyPI version](https://badge.fury.io/py/Truvari.svg)](https://badge.fury.io/py/Truvari)
[![pylint](imgs/pylint.svg)](https://github.com/acenglish/truvari/actions/workflows/pylint.yml)
[![FuncTests](https://github.com/acenglish/truvari/actions/workflows/func_tests.yml/badge.svg?branch=develop&event=push)](https://github.com/acenglish/truvari/actions/workflows/func_tests.yml)
[![coverage](imgs/coverage.svg)](https://github.com/acenglish/truvari/actions/workflows/func_tests.yml)
[![develop](https://img.shields.io/github/commits-since/acenglish/truvari/v5.2.0)](https://github.com/ACEnglish/truvari/compare/v5.2.0...develop)
[![Downloads](https://static.pepy.tech/badge/truvari)](https://pepy.tech/project/truvari)

![Logo](https://raw.githubusercontent.com/ACEnglish/truvari/develop/imgs/BoxScale1_DarkBG.png)  
Toolkit for benchmarking, merging, and annotating Structural Variants

📚 [WIKI page](https://github.com/acenglish/truvari/wiki) has detailed user documentation.  
🛠️ [Developer Docs](https://truvari.readthedocs.io/en/latest/) for the truvari API.  
📈 See [Updates](https://github.com/acenglish/truvari/wiki/Updates) on new versions.  
📝 Read our Papers ([#1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02840-6), [#2](https://rdcu.be/dFQNN)) to learn more.

## 💻 Installation
Truvari uses Python 3.6+ and can be installed with pip:
```
  python3 -m pip install Truvari 
```
For details and more installation options, see [Installation](https://github.com/acenglish/truvari/wiki/Installation) on the wiki.

## ⏩ Quick Start

Each sub-command contains help documentation. Start with `truvari -h` to see available commands.

The current most common Truvari use case is for structural variation benchmarking:
```
  truvari bench -b base.vcf.gz -c comp.vcf.gz -f reference.fa -o output_dir/
```

Find more matches by harmonizing phased variants using refine:
```
   truvari refine output_dir/
```

Use Truvari's comparison engine to consolidate redundant variants in a merged multi-sample VCF:
```
    bcftools merge -m none sampleA.vcf.gz sampleB.vcf.gz | bgzip > merge.vcf.gz
    tabix merge.vcf.gz
    truvari collapse -i merge.vcf.gz -o truvari_merge.vcf
```

### truvari survivor
N-way VCF merging with SURVIVOR-style support tracking. This command groups variants from multiple VCFs and annotates the merged records with support metrics.

```
    truvari survivor -i caller1.vcf caller2.vcf caller3.vcf -o merged.vcf
```

**Key Parameters:**
- `-i, --input`: List of input VCFs to merge.
- `-l, --list`: A text file containing a list of VCF paths to merge.
- `-o, --output`: Path to the output merged VCF.
- `-r, --refdist`: Max reference location distance for grouping (default: 500).
- `-P, --pctsize`: Min percent size similarity (default: 0.50).
- `-O, --pctovl`: Min reciprocal overlap for non-BND variants (default: 0.50).
- `-p, --pctseq`: Min percent sequence similarity (default: 0.0).
- `-B, --bnddist`: Max distance allowed between BNDs (defaults to refdist).
- `--strandignore`: Ignore strand/orientation mismatches during BND comparison (default: True).
- `-D, --decompose`: Allow decomposition of symbolic SVs to BNDs for cross-type comparison.
- `-s, --sizemin`: Minimum variant size to consider (default: 50).
- `--passonly`: Only consider variants with `FILTER == PASS`.

**Merged Annotations:**
- `SUPP`: Number of callers supporting the variant.
- `SUPP_VEC`: A bit-vector (string) indicating which callers support the variant.
- `CALLERS`: Names of the supporting callers.
- `CALLER_IDS`: Original variant IDs from each supporting caller.

## 🧬 Truvari Commands

 - [bench](https://github.com/acenglish/truvari/wiki/bench) - Performance metrics from comparison of two VCFs
 - [survivor](https://github.com/acenglish/truvari/wiki/survivor) - N-way VCF merging with support tracking
 - [collapse](https://github.com/acenglish/truvari/wiki/collapse) - Collapse possibly redundant VCF entries
 - [refine](https://github.com/ACEnglish/truvari/wiki/refine) - Automated bench result refinement with phab
 - [anno](https://github.com/acenglish/truvari/wiki/anno) - Add SV annotations to a VCF
 - [phab](https://github.com/ACEnglish/truvari/wiki/phab) - Harmonize variant representations using MSA
 - [consistency](https://github.com/acenglish/truvari/wiki/consistency) - Consistency report between multiple VCFs
 - [vcf2df](https://github.com/acenglish/truvari/wiki/vcf2df) - Turn a VCF into a pandas DataFrame
 - [segment](https://github.com/acenglish/truvari/wiki/segment) - Normalization of SVs into disjointed genomic regions
 - [stratify](https://github.com/acenglish/truvari/wiki/stratify) - Count variants per-region in vcf
 - [divide](https://github.com/ACEnglish/truvari/wiki/divide) - Divide a VCF into independent shards
 - [ga4gh](https://github.com/ACEnglish/truvari/wiki/ga4gh) - Consolidate benchmarking result VCFs

## 🔎 More Information

All documentation about Truvari is on the [WIKI](https://github.com/acenglish/truvari/wiki). Additional information about using Truvari can be found in [Discussions](https://github.com/ACEnglish/truvari/discussions)
