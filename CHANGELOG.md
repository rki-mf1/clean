# Changelog

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
