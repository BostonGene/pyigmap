# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

* Add a new flag `--discard-junctions-with-N` toto `cdr3nt-error-corrector` step ([#92](https://github.com/BostonGene/pyigmap/issues/92))
* Podman support ([#69](https://github.com/BostonGene/pyigmap/issues/41))
* Add a new flag `--only-canonical` to `cdr3nt-error-corrector` step ([#41](https://github.com/BostonGene/pyigmap/issues/41)) ([217c316](https://github.com/BostonGene/pyigmap/commit/217c316f82a9613a0b3e5994f90b50fbed3e37b6))
* Add `--filter-pgen-singleton`, `--filter-pgen-all` and `--skip-pgen-calculation` flags to `cdr3nt-error-corrector` ([#48](https://github.com/BostonGene/pyigmap/issues/48)) ([f04b3cd](https://github.com/BostonGene/pyigmap/commit/f04b3cd36646e0a7272d3026794c4b726f59d7af))
* Automate python, java and docker installation ([#53](https://github.com/BostonGene/pyigmap/issues/53)) ([25c8bf8](https://github.com/BostonGene/pyigmap/commit/25c8bf8e0dd68a872acd77f1648d2bc92ad15ec7))
* New metrics: `no_v_call` and `no_j_call` in `stat.json` (`cdr3nt-error-corrector` step) ([#40](https://github.com/BostonGene/pyigmap/issues/40)) ([c064525](https://github.com/BostonGene/pyigmap/commit/c0645254485e62a26f6d17aa913d632bc751a8d7))
* Use `nextflow` to build pipeline ([#2](https://github.com/BostonGene/pyigmap/issues/2)) ([b575887](https://github.com/BostonGene/pyigmap/commit/b57588781ea7cd4aaaaac869d6fe4df041159da1))
* Support `pre-commit` hooks and `Makefile` shortcuts ([#27](https://github.com/BostonGene/pyigmap/issues/27)) ([c1ec4be](https://github.com/BostonGene/pyigmap/commit/c1ec4be848a875335f6aaa03535eeeca9ee734a1))
* Support downloading test datasets from [zenodo](https://zenodo.org/records/11103555) by `--zenodo` flag ([#38](https://github.com/BostonGene/pyigmap/issues/38)) ([77bdf8c](https://github.com/BostonGene/pyigmap/commit/77bdf8c5381370a40daabaa09ff1760b1ea36770))

### Changed

* Consider mock merging reads for non-overlapped reads ([#84](https://github.com/BostonGene/pyigmap/issues/84))
* Pgen calculation disabled for amplicon ([#76](https://github.com/BostonGene/pyigmap/issues/76))
* Vidjil disabled for amplicon data ([#75](https://github.com/BostonGene/pyigmap/issues/75))
* Optimized workflow: igblast processes TCR and BCR in one run ([#78](https://github.com/BostonGene/pyigmap/issues/78))
* Now TCR and BCR run in one vidjil job ([#67](https://github.com/BostonGene/pyigmap/issues/67))
* Remove clones with `junction == None` and disable `_process_cdr3_sequences()` function for annotation generated using IgBLAST ([#80](https://github.com/BostonGene/pyigmap/issues/80))
* Remove the same filters in cdr3nt-error-corrector ([#60](https://github.com/BostonGene/pyigmap/issues/60)) ([fff8937](https://github.com/BostonGene/pyigmap/commit/fff8937f6e8c420b4cb2c3409dcd13c530e056ca))
* Reference generation via Makefile ([#55](https://github.com/BostonGene/pyigmap/issues/55)) ([03824c](https://github.com/BostonGene/pyigmap/commit/803824ca3479fd121802281bea4071cd719230c0))
* Move most of the functions from `run.py` to `utils.py` (for compatibility with internal bitbucket repository) ([#28](https://github.com/BostonGene/pyigmap/issues/28)) ([d3fdfd7](https://github.com/BostonGene/pyigmap/commit/d3fdfd71e33f2995bb349c7880bd9d5af06c4d7a))
* Make Dockerfiles [multi-staged](https://docs.docker.com/build/building/multi-stage/) to build an executable and not executable images ([#44](https://github.com/BostonGene/pyigmap/issues/44)) ([c602688](https://github.com/BostonGene/pyigmap/commit/c602688ac17bc1259f52394eebaac7e83167d459))
* Select *01 (major) allele as a default in the igblast ref ([#26](https://github.com/BostonGene/pyigmap/issues/26)) ([db50abe](https://github.com/BostonGene/pyigmap/commit/db50abedf31e47d2b9f0791e6c653dd9a6e0f732))

### Fixed

* Fix a bug with chimeras removing ([#49](https://github.com/BostonGene/pyigmap/issues/49)) ([f80356d](https://github.com/BostonGene/pyigmap/commit/f80356de9d3ed37dad42d2c554f76c488050760c))
* Fix a bug with dataframe index in `drop_duplicates_in_different_loci` ([#29](https://github.com/BostonGene/pyigmap/issues/29)) ([dda232c](https://github.com/BostonGene/pyigmap/commit/dda232c20c28c72403cc3a08db90d6fe33620c85))
* Fix a bug with `total_reads` calculation for fastp.json in `cdr3nt-error-corrector` step ([#16](https://github.com/BostonGene/pyigmap/issues/16)) ([16fb28b](https://github.com/BostonGene/pyigmap/commit/16fb28b81477545fb392d96b05d431b85bd2a0d4))
* Fix an `igblast` reference ([#11](https://github.com/BostonGene/pyigmap/pull/11)) ([d3511d5](https://github.com/BostonGene/pyigmap/commit/d3511d57d9856f88a49e4c884b9d1650bc091d18))
