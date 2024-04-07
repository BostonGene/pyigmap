install: ## Install and check dependencies
	@java --version
	@docker version
	curl -s https://get.nextflow.io | bash
	docker build -t calib-dedup steps/calib_dedup
	docker build -t fastp steps/fastp
	docker build -t vidjil steps/vidjil
	docker build -t igblast steps/igblast
	docker build -t cdr3nt-error-corrector steps/cdr3nt_error_corrector
	docker build -t downloader steps/downloader
