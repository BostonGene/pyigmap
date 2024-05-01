# Changelog

## Unreleased

### Features

* Use `nextflow` to build pipeline ([#2](https://github.com/BostonGene/pyigmap/issues/2))
* Add new flag `--calculate-pgen` to `igblast` step [#3](https://github.com/BostonGene/pyigmap/issues/3) ([eedf2b2](https://github.com/BostonGene/pyigmap/commit/eedf2b20487bbb68882b4e56631f72c664ec4167))
* Support `pre-commit` hooks and `Makefile` shortcuts ([#27](https://github.com/BostonGene/pyigmap/issues/27)) ([c1ec4be](https://github.com/BostonGene/pyigmap/commit/c1ec4be848a875335f6aaa03535eeeca9ee734a1))
* Add new `--filter-singleton` flag to `cdr3nt-error-corrector step` ([#30](https://github.com/BostonGene/pyigmap/issues/30)) ([40acf9e](https://github.com/BostonGene/pyigmap/commit/40acf9ee0ec29e0a6e901123c2be2b0c99a70654))

### Bug Fixes

* Fix bug with dataframe index in `drop_duplicates_in_different_loci` ([#29](https://github.com/BostonGene/pyigmap/issues/29)) ([dda232c](https://github.com/BostonGene/pyigmap/commit/dda232c20c28c72403cc3a08db90d6fe33620c85))
* Fix bug with `total_reads` calculation for fastp.json in `cdr3nt-error-corrector` step ([#16](https://github.com/BostonGene/pyigmap/issues/16)) ([16fb28b](https://github.com/BostonGene/pyigmap/commit/16fb28b81477545fb392d96b05d431b85bd2a0d4))
* Fix `igblast` reference ([#11](https://github.com/BostonGene/pyigmap/pull/11)) ([d3511d5](https://github.com/BostonGene/pyigmap/commit/d3511d57d9856f88a49e4c884b9d1650bc091d18))