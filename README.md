# SwiftCov

SwiftCov is a drop-in replacement for [bedtools](https://github.com/arq5x/bedtools2) v2.31.1 `coverage`. It writes
the same output, byte for byte, much faster.

On real HG002 reads with every query row on the chromosome that has reads, it was 13- to 34-fold faster than bedtools
in its faster mode, `-sorted`, at one thread (query restricted to chr21: 13-fold on exome targets, 34-fold on
transcripts). It was 9.2-fold faster than bedtk, the fastest other tool, at one thread (6- to 11-fold on four tasks,
genome-wide queries), and bedtk writes a different format. No output byte differed from bedtools in 22,200 generated
test cases or 4,194 benchmark runs. Details are in the preprint (see Citation).

The larger factors are for genome-wide queries, where most rows lie on chromosomes with no reads, which inflates them:

| Genome-wide query, against bedtools `-sorted` | 1 thread | 2 threads |
|---|---:|---:|
| 12.5 million chr21 alignments | 36-fold | 50-fold |
| 10.5 million chr22 alignments, held out | 40-fold | 49-fold |

These are geometric means on a 2-vCPU cloud VM. On an AVX2 laptop at 12 threads the genome-wide figures were 25- and
38-fold.

## Build

Requirements: Linux on x86-64, git, make, gcc/g++ with C++17 (tested with 13.3.0 on Ubuntu 24.04), Python 3 (standard
library only), and network access to GitHub on the first build, which clones zlib v1.3.1 into `ext/deps` and builds it
as a static library.

```
sudo apt install build-essential git python3
bash solver/build.sh        # produces solver/swiftcov; SWIFTCOV_BUILD_JOBS=N sets make -j (default 2)
```

AVX2 is used when the CPU has it; other x86-64 CPUs use SSE2 code paths and give the same output, more slowly.

## Usage

`solver/swiftcov` is a complete bedtools v2.31.1 binary. Use its `coverage` subcommand exactly as `bedtools coverage`,
with the same flags:

```
solver/swiftcov coverage -a targets.bed -b reads.bed > result.txt
OMP_NUM_THREADS=2 solver/swiftcov coverage -a targets.bed -b reads.bed -hist -sorted > result.txt
```

- Standard output and created files are byte-identical to bedtools v2.31.1, at every thread count, and where bedtools
  exits non-zero SwiftCov does too. Stderr wording is not compared.
- `OMP_NUM_THREADS` sets the number of threads (default 1).
- Fast path, where the speedup comes from: `-a` with one file and `-b` with one or more files, all tab-delimited BED
  (BED3 to BED6 and BED-plus with a constant field count) with canonical integer coordinates and no header,
  carriage-return, NUL or empty lines; all five output modes (default, `-counts`, `-mean`, `-hist`, `-d`); `-sorted` and
  `-g`; `-s`/`-S` on BED6 or wider; `-f`, `-F`, `-r`, `-e` (except `-e` with `-f` and no `-F`); and the output-neutral
  `-header`, `-nonamecheck`, `-nobuf`, `-bed`, `-sortout`, `-ubam`, `-fbam`, `-seed`.
- Everything else runs the unmodified bedtools code in a fresh process, with the same output and no speedup: gzip or
  BGZF input, BAM, GFF/GTF, VCF, BED12 with `-split`, files with headers or CRLF line endings, non-canonical numbers,
  other options and unsorted `-sorted` input. Other subcommands are also the unmodified bedtools code.
- Not tested: CRAM, remote URLs and reading from standard input.
- Memory: SwiftCov maps the read file and peaked at 548 to 625 MiB on these sets, against 8 to 332 MiB for bedtools
  `-sorted` and 3.8 to 4.5 GiB for its default mode.
- `SWIFTCOV_DEBUG=1` prints why an input went to the bedtools code; `SWIFTCOV_NOFAST=1` turns the fast path off.

## Check it yourself

```
python3 prepare.py              # builds bedtools v2.31.1 (commit 705ccfd), rebuilds the benchmark sets from data/ and
                                # records bedtools' outputs and times (needs gcc, g++, git, network)
python3 evaluate.py             # builds SwiftCov; 300 fresh test cases and three benchmark sets compared byte for byte
python3 evaluate.py --holdout   # the same on the four holdout sets
```

`prepare.py` needs about 6 GB of free memory, because bedtools' default mode holds 4.5 GiB of reads, and about 3 GB of
disk. It took about 10 minutes on the 2-vCPU VM. `evaluate.py` prints `gate: PASS` and the speedups for genome-wide
queries; add `--threads 1` for one thread. `real_speedup_vs_upstream_best` is the comparison in the table above.
Timings depend on the host.

For the query restricted to chr21, after the commands above:

```
awk '$1 == "chr21"' bench/data/source_refseq_transcripts.bed > chr21.bed
time ext/bedtools/bin/bedtools coverage -a chr21.bed -b bench/data/source_dev_reads.bed -sorted > bedtools.txt
time OMP_NUM_THREADS=1 solver/swiftcov coverage -a chr21.bed -b bench/data/source_dev_reads.bed > swiftcov.txt
cmp bedtools.txt swiftcov.txt
```

Use `source_exome_targets.bed` and add `-d` to both tools for the exome case. `data/` holds the real reads and
annotations used in the preprint as compact lossless encodings (see `data/README.txt`).

## Licence

MIT (`LICENSE`), the licence of bedtools v2.31.1. bedtools' `LICENSE` moved to MIT in 2019; 105 of its source files
still carry an older GPL 2.0 header, which are left as they are; `LICENSE` states what applies. Vendored with bedtools:
htslib (MIT/Expat, modified BSD for its `cram/` directory) and gzstream (LGPL-2.1). The build links zlib statically (zlib
licence). The data in `data/` are not covered by `LICENSE`: the HG002 reads are CC0, and the RefSeq and Twist annotations
are redistributed under UCSC's data-reuse terms and NCBI's policy for public molecular data (`data/README.txt`).

## Citation

Christensen J. SwiftCov computes read coverage 13- to 34-fold faster than bedtools, with byte-identical output.
Preprint, 2026. Please also cite bedtools (Quinlan AR, Hall IM, Bioinformatics 26:841–842,
2010). See `CITATION.cff`.
