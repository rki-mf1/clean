# Changelog

## unreleased

### Changed

- replaced `bwa mem` with [`bwa-mem2`](https://github.com/bwa-mem2/bwa-mem2) as short-read mapper alternative; the `--bwa` parameter is unchanged
  - faster index building, which mainly helps for large indices combining several eukaryotic genomes
  - note that a `bwa-mem2` index is about 4x larger than a `bwa` index and is held in memory during mapping, so the memory requirements went up

- all containers and conda environments now ship the same `samtools`/`htslib` version (1.24); this also updates `minimap2` 2.26 -> 2.31, `bedtools` 2.30.0 -> 2.31.1, `seqkit` 2.6.1 -> 2.13.0 and `pigz` 2.3.4 -> 2.8

### Fixed

- the short-read mapper alternative now also gets an environment with the `conda`/`mamba` profiles and CPUs/memory with the `local`/`standard` profiles
- `samclipy` gets its own container (the `samtools` one, it has `python` and `git`) instead of implicitly using the one of the `smallTask` label, and `git` was added to its conda environment
- syntax that the strict parser of Nextflow >=25.10 rejects: the variable declaration in `nextflow.config`, typed `for` loops, `if` blocks around `publishDir`/`storeDir` directives, `env(VAR)` outputs and `addParams()` on `include` statements
  - the CI matrix is pinned to Nextflow 25.10.0 instead of `latest`, because the top-level statements in `clean.nf` still need to move into the entry workflow before the pipeline runs on Nextflow >=26.04

## [v1.1.0] - 2024-11-08

### Added

- `bwa mem` as short-read mapper alternative, parameter: `--bwa`

## [v1.0.3] - 2024-08-08

### Added

- Set default branch to `main` instead of the Nextflow default `master`
- Remove `conda clean` from GitHub action to avoid random crashes
- Bump github action versions for node 16 -> 20 change
- Add T2T homo sapiens genome as additional auto-download option

## [v1.0.2] - 2024-05-17

### Added

- added `--skip_qc` option to skip QC steps

## [v1.0.1] - 2024-03-15

### Added

- SARS-CoV-2 added to the auto-download option

## [v1.0.0] - 2024-01-04

### Changed

- reorganized of results directory and file names

### Added

- dry run CI tests
- options to reduce disk usage
  - `--cleanup_work_dir` and `--no_intermediate`

### Fixed

- fixed some issues on Mac OS

## [v1.0.0-beta.1] - 2023-10-11

### Changed

- changed minimap2 container so that `ncurses` is included

## [v1.0.0-beta] - 2023-09-30

### Changed

- changed input parameter usage:
  - before: `--[nano|illumina|illumina_single_end|fasta]` 
  - now: `--input_type [nano|illumina|illumina_single_end|fasta] --input *.fastq`
- changed workflow figure to a nicer figure
- changed workflow structure (introducing subworkflows)
- input files with the suffix `clean` are not allowed 

### Added

- added CHANGELOG.md, Citations.md and citation information
- added `--cleanup_work_dir` to remove work dir files after a successful run
- added `--min_clip` to filter mapped reads by soft-clipped length
- added `--dcs_strict` to use only DCS reads with artificial ends
- added `stub` command for Nextflow prototyping
- added `idxstats` 

## Fixed

- pipeline report with timestamp
- `--split-prefix` parameter for `minimap2`
- make concat contamination more efficient
