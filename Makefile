SHELL:=/bin/bash
VIRTUAL_ENV=env
PYTHON=${VIRTUAL_ENV}/bin/python3
STAGE=not_exec

# .ONESHELL:
DEFAULT_GOAL: install
.PHONY: help run clean build test test_wf mypy check format update venv

# Colors for echos 
ccend = $(shell tput sgr0)
ccbold = $(shell tput bold)
ccgreen = $(shell tput setaf 2)
ccso = $(shell tput smso)

mypy: venv ## >> run mypy type checker
	@echo ""
	@echo "$(ccso)--> Running mypy $(ccend)"
	$(PYTHON) -m mypy steps/calib_dedup/
	$(PYTHON) -m mypy steps/fastp/
	$(PYTHON) -m mypy steps/vidjil/
	$(PYTHON) -m mypy steps/igblast/
	$(PYTHON) -m mypy steps/cdr3nt_error_corrector/

check: venv ## >> run ruff linter
	@echo ""
	@echo "$(ccso)--> Running ruff check $(ccend)"
	$(PYTHON) -m ruff check steps/calib_dedup/ steps/fastp/ steps/vidjil/ steps/igblast/ steps/cdr3nt_error_corrector/

format: venv ## >> run ruff formatter
	@echo ""
	@echo "$(ccso)--> Running ruff format $(ccend)"
	$(PYTHON) -m ruff format steps/calib_dedup/ steps/fastp/ steps/vidjil/ steps/igblast/ steps/cdr3nt_error_corrector/

test_wf: venv ## >> run tests for all workflows via pytest and pytest-workflow tool
	@echo ""
	@echo "$(ccso)--> Running workflow tests $(ccend)"
	$(PYTHON) -m pytest tests/ -vv

test: venv ## >> run tests for all steps via pytest tool
	@echo ""
	@echo "$(ccso)--> Running steps tests $(ccend)"
#	$(PYTHON) -m pytest steps/calib_dedup/unit_tests -vv
#	$(PYTHON) -m pytest steps/fastp/unit_tests -vv
#	$(PYTHON) -m pytest steps/vidjil/unit_tests -vv
	$(PYTHON) -m pytest steps/igblast/unit_tests -vv
	$(PYTHON) -m pytest steps/cdr3nt_error_corrector/unit_tests -vv

clean: ## >> remove docker images, python environment and nextflow build files
	@echo ""
	@echo "$(ccso)--> Removing virtual environment $(ccend)"
	docker rmi -f downloader calib-dedup fastp vidjil igblast cdr3nt-error-corrector
	rm -rf $(VIRTUAL_ENV) .nextflow.log* work .nextflow nextflow

build: ##@main >> build docker images, the virtual environment and install requirements
	@echo ""
	@echo "$(ccso)--> Build $(ccend)"
	$(MAKE) clean
	$(MAKE) STAGE=exec
	$(MAKE) update

venv: $(VIRTUAL_ENV)

$(VIRTUAL_ENV): ## >> install virtualenv and setup the virtual environment
	@echo "$(ccso)--> Install and setup virtualenv $(ccend)"
	python3 -m pip install --upgrade pip
	python3 -m pip install virtualenv
	virtualenv $(VIRTUAL_ENV)

update: venv ##@main >> update requirements.txt inside the virtual environment
	@echo "$(ccso)--> Updating packages $(ccend)"
	$(PYTHON) -m pip install -r ./steps/calib_dedup/requirements.txt
	$(PYTHON) -m pip install -r ./steps/igblast/requirements.txt
	$(PYTHON) -m pip install -r ./steps/cdr3nt_error_corrector/requirements.txt
	$(PYTHON) -m pip install pytest==8.1.1 pytest-workflow==2.1.0 ruff==0.4.2 mypy==1.10.0

install: ## Install and check dependencies
	@java --version
	@docker version
	curl -s https://get.nextflow.io | bash
	docker build --target $(STAGE) -t calib-dedup steps/calib_dedup
	docker build --target $(STAGE) -t fastp steps/fastp
	docker build --target $(STAGE) -t vidjil steps/vidjil
	docker build --target $(STAGE) -t igblast steps/igblast
	docker build --target $(STAGE) -t cdr3nt-error-corrector steps/cdr3nt_error_corrector
	docker build -t downloader steps/downloader

# And add help text after each target name starting with '\#\#'
# A category can be added with @category
HELP_FUN = \
	%help; \
	while(<>) { push @{$$help{$$2 // 'options'}}, [$$1, $$3] if /^([a-zA-Z\-\$\(]+)\s*:.*\#\#(?:@([a-zA-Z\-\)]+))?\s(.*)$$/ }; \
	print "usage: make [target]\n\n"; \
	for (sort keys %help) { \
	print "${WHITE}$$_:${RESET}\n"; \
	for (@{$$help{$$_}}) { \
	$$sep = " " x (32 - length $$_->[0]); \
	print "  ${YELLOW}$$_->[0]${RESET}$$sep${GREEN}$$_->[1]${RESET}\n"; \
	}; \
	print "\n"; }

help: ##@other >> Show this help.
	@perl -e '$(HELP_FUN)' $(MAKEFILE_LIST)
	@echo ""
	@echo "Note: to activate the environment in your local shell type:"
	@echo "   $$ source $(VIRTUAL_ENV)/bin/activate"
