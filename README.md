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
2. Run `gffquant <input_gff> <input_bam> -o <out_prefix> [--ambig_mode {[unique_only], all1, 1overN}]
  - The `<input_bam>` file needs to be position sorted.
  - Output files are `<out_prefix>.seqname.txt` (contig counts) and `<out_prefix>.feature_counts.txt` (feature/subfeature counts).
  - `--ambig_mode` controls how ambiguously mapping reads are processed. These are analogous to the ngless modes:
      - `'unique_only` will simply ignore any reads that is labeled as ambiguous (`MAPQ=0`). Runtime increases with bamfile size, memory should remain more or less constant.
	  - `all1` will treat all ambiguous alignments as single ended individual reads. Runtime increases with bamfile size, memory should remain more or less constant.
	  - `1overN` will dump all ambiguous alignments to disk (in reduced form), finish the processing of the unique alignments, then read and process the ambiguous alignments. Runtime and memory increase significantly with bamfile size. In addition, temporary disk space proportional to the number of ambiguous alignments is needed. In this mode (and only in this mode), sequence counting will include an additional output with reads distributed to reference contigs according to the `dist1` mode of ngless.
	  
### Resource requirements
Memory and computing time requirements correspond to bamfile size. For current bamfile sizes of up to 24GB:
  - `uniq_only` and `all1`: ~10GB memory and ~4h
  - `1overN`: >10GB memory, ~10min - >8h
  
### TODO
