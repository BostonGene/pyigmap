# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

* New metrics: `no_v_call` and `no_j_call` in `stat.json` (`cdr3nt-error-corrector` step) ([#40](https://github.com/BostonGene/pyigmap/issues/40)) ([c064525](https://github.com/BostonGene/pyigmap/commit/c0645254485e62a26f6d17aa913d632bc751a8d7))
* Use `nextflow` to build pipeline ([#2](https://github.com/BostonGene/pyigmap/issues/2)) ([b575887](https://github.com/BostonGene/pyigmap/commit/b57588781ea7cd4aaaaac869d6fe4df041159da1))
* Add new flag `--calculate-pgen` to `igblast` step [#3](https://github.com/BostonGene/pyigmap/issues/3) ([eedf2b2](https://github.com/BostonGene/pyigmap/commit/eedf2b20487bbb68882b4e56631f72c664ec4167))
* Support `pre-commit` hooks and `Makefile` shortcuts ([#27](https://github.com/BostonGene/pyigmap/issues/27)) ([c1ec4be](https://github.com/BostonGene/pyigmap/commit/c1ec4be848a875335f6aaa03535eeeca9ee734a1))
* Add new `--filter-singleton` flag to `cdr3nt-error-corrector` step ([#30](https://github.com/BostonGene/pyigmap/issues/30)) ([40acf9e](https://github.com/BostonGene/pyigmap/commit/40acf9ee0ec29e0a6e901123c2be2b0c99a70654))
* Support downloading test datasets from [zenodo](https://zenodo.org/records/11103555) by `--zenodo` flag ([#38](https://github.com/BostonGene/pyigmap/issues/38)) ([77bdf8c](https://github.com/BostonGene/pyigmap/commit/77bdf8c5381370a40daabaa09ff1760b1ea36770))

### Changed

* Move most of the functions from `run.py` to `utils.py` (for compatibility with internal bitbucket repository) ([#28](https://github.com/BostonGene/pyigmap/issues/28)) ([d3fdfd7](https://github.com/BostonGene/pyigmap/commit/d3fdfd71e33f2995bb349c7880bd9d5af06c4d7a))
* Make Dockerfiles [multi-staged](https://docs.docker.com/build/building/multi-stage/) to build an executable and not executable images ([#44](https://github.com/BostonGene/pyigmap/issues/44)) ([c602688](https://github.com/BostonGene/pyigmap/commit/c602688ac17bc1259f52394eebaac7e83167d459))
* Select *01 (major) allele as a default in the igblast ref ([#26](https://github.com/BostonGene/pyigmap/issues/26)) ([db50abe](https://github.com/BostonGene/pyigmap/commit/db50abedf31e47d2b9f0791e6c653dd9a6e0f732))

### Fixed

* Fix bug with dataframe index in `drop_duplicates_in_different_loci` ([#29](https://github.com/BostonGene/pyigmap/issues/29)) ([dda232c](https://github.com/BostonGene/pyigmap/commit/dda232c20c28c72403cc3a08db90d6fe33620c85))
* Fix bug with `total_reads` calculation for fastp.json in `cdr3nt-error-corrector` step ([#16](https://github.com/BostonGene/pyigmap/issues/16)) ([16fb28b](https://github.com/BostonGene/pyigmap/commit/16fb28b81477545fb392d96b05d431b85bd2a0d4))
* Fix `igblast` reference ([#11](https://github.com/BostonGene/pyigmap/pull/11)) ([d3511d5](https://github.com/BostonGene/pyigmap/commit/d3511d57d9856f88a49e4c884b9d1650bc091d18))
