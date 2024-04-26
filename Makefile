SHELL:=/bin/bash
VIRTUAL_ENV=env
PYTHON=${VIRTUAL_ENV}/bin/python3

# .ONESHELL:
DEFAULT_GOAL: help
.PHONY: help run clean build test

# Colors for echos 
ccend = $(shell tput sgr0)
ccbold = $(shell tput bold)
ccgreen = $(shell tput setaf 2)
ccso = $(shell tput smso)

test: ## >> running tests for all steps via pytest tool
	@echo ""
	@echo "$(ccso)--> Running tests $(ccend)"
#	. $(VIRTUAL_ENV)/bin/activate; python3 -m pytest steps/calib_dedup/unit_tests -v
#	. $(VIRTUAL_ENV)/bin/activate; python3 -m pytest steps/fastp/unit_tests -v
#	. $(VIRTUAL_ENV)/bin/activate; python3 -m pytest steps/vidjil/unit_tests -v
#	. $(VIRTUAL_ENV)/bin/activate; python3 -m pytest steps/igblast/unit_tests -v
	. $(VIRTUAL_ENV)/bin/activate; python3 -m pytest steps/cdr3nt_error_corrector/unit_tests -v

clean: ## >> remove all environment and build files
	@echo ""
	@echo "$(ccso)--> Removing virtual environment $(ccend)"
	rm -rf $(VIRTUAL_ENV)

build: ##@main >> build the virtual environment with an ipykernel for jupyter and install requirements
	@echo ""
	@echo "$(ccso)--> Build $(ccend)"
	$(MAKE) clean
	$(MAKE) install

venv: $(VIRTUAL_ENV) ## >> install virtualenv and setup the virtual environment

$(VIRTUAL_ENV):
	@echo "$(ccso)--> Install and setup virtualenv $(ccend)"
	python3 -m pip install --upgrade pip
	python3 -m pip install virtualenv
	virtualenv $(VIRTUAL_ENV)

install: venv ##@main >> update requirements.txt inside the virtual environment
	@echo "$(ccso)--> Updating packages $(ccend)"
	$(PYTHON) -m pip install -r ./steps/calib_dedup/requirements.txt
	$(PYTHON) -m pip install -r ./steps/igblast/requirements.txt
	$(PYTHON) -m pip install -r ./steps/cdr3nt_error_corrector/requirements.txt
	$(PYTHON) -m pip install pytest==8.1.1


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