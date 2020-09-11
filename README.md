### Prerequisites
- Python3.7+ (might work with Python3.6, but not tested)

### Dependencies
(will be installed via pip)
- pyyaml
- intervaltree

### Installation
1. Obtain the source code `git clone https://github.com/cschu/gff_quantifier.git`.
2. `cd gff_quantifier`.
3. `python setup.py bdist_wheel`.
4. `pip install [--user] -U dist/gffquant-<version>-py3-none-any.whl`.

### Running gffquant
1. Pre-index your gff-file with `gff_indexer <input_gff>`.
  - This will write the index to `<input_gff>.index` and only needs to be done once per gff.
  - For best results, the gff should be strictly sorted by seqname (column 1).
  - The gff must not be gzipped (as it is otherwise not feasible to seek()).
2. Run `gffquant <input_gff> <input_bam> -o <out_prefix>` to compute the read counts against the provided gff.
  - `<input_gff>.index` needs to be in the same path as `<input_gff>`.
  - The bam file needs to be position sorted.
  - Memory and computing time requirements correspond to bamfile size (current bamfile sizes of up to 24GB should be fine with 8GB memory and 4h.)
  - This will generate the files `<out_prefix>.seqname.txt` (contig counts) and `<out_prefix>.feature_counts.txt` (feature/subfeature counts).
  
### Caveats
- The current version (0.2) only quantifies unique, non-supplementary alignments (support for ambiguous alignments is being developed, 
but will require an additional name-sorted bamfile in order to keep memory requirements for large samples down.)
