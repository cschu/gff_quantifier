### Prerequisites
- Python3.7+ (we need to be able to rely on dictionary item order preservation!)

### Dependencies
(will be installed via pip)
- pyyaml
- intervaltree

### Installation
1. Obtain the source code `git clone https://github.com/cschu/gff_quantifier.git`.
2. `cd gff_quantifier`.
3. `python setup.py bdist_wheel`.
4. `pip install [--user] -U dist/gffquant-<version>-py3-none-any.whl`.

After this, the relevant commands `gff_indexer` and `gffquant` should be in your path.

### Running gffquant
1. Pre-index your gff-file with `gff_indexer <input_gff>`.
  - This will write the index to `<input_gff>.index` and only needs to be done once per gff.
  - For best results, the gff should be strictly sorted by seqname (column 1).
  - The gff must not be gzipped (random access of gzipped files via seek() is not feasible, hence gzipped gffs are not supported).
2. (optional) Produce a reduced bamfile, containing only the ambiguous alignments.
  - `samtools view -h -F 0x800 <input_bam> | awk '/^[^@]/ { if ($5 != 0) next; } { print $0; } | samtools sort -@ <threads> -n -o <name_sorted_bam>`
3. Run `gffquant <input_gff> <input_bam> -o <out_prefix> [-n <name_sorted_bam>]` to compute the read counts against the provided gff.
  - `<input_gff>.index` needs to be in the same path as `<input_gff>`.
  - The `<input_bam>` file needs to be position sorted, the `<name_sorted_bam>` file by name.
  - Memory and computing time requirements correspond to bamfile size. For current bamfile sizes of up to 24GB:
    - unique alignments only: ~10GB memory and ~4h
	- including ambiguous reads: ~10GB memory and >10h (not set in stone yet)
  - This will generate the files `<out_prefix>.seqname.txt` (contig counts) and `<out_prefix>.feature_counts.txt` (feature/subfeature counts).
  
### Caveats
- Processing of ambiguous alignments is excruciatingly slow in the current version. I am looking into distributing the bam processing onto multiple threads.
