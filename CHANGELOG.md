# Changelog

## unreleased

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
