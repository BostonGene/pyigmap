SHELL:=/bin/bash

PYTHON_VERSION=3.9.19
PYTHON_SYS=python3.9
VIRTUAL_ENV=env
PYTHON_ENV=${VIRTUAL_ENV}/bin/python3
JAVA_VERSION=21
ARCHITECTURE=amd64
BUILD_REF_STAGE=build-ref

# .ONESHELL:
DEFAULT_GOAL: install
.PHONY: help clean build build-ref integration-tests unit-tests mypy check format install-nextflow install-python install-docker install-java

# Colors for echos 
ccend = $(shell tput sgr0)
ccbold = $(shell tput bold)
ccgreen = $(shell tput setaf 2)
ccso = $(shell tput smso)

mypy: venv ## >> run mypy type checker
	@echo ""
	@echo "$(ccso)--> Running mypy $(ccend)"
	$(PYTHON_ENV) -m mypy steps/calib_dedup/
	$(PYTHON_ENV) -m mypy steps/fastp/
	$(PYTHON_ENV) -m mypy steps/vidjil/
	$(PYTHON_ENV) -m mypy steps/igblast/
	$(PYTHON_ENV) -m mypy steps/cdr3nt_error_corrector/

check: venv ## >> run ruff linter
	@echo ""
	@echo "$(ccso)--> Running ruff check $(ccend)"
	$(PYTHON_ENV) -m ruff check steps/calib_dedup/ steps/fastp/ steps/vidjil/ steps/igblast/ steps/cdr3nt_error_corrector/

format: venv ## >> run ruff formatter
	@echo ""
	@echo "$(ccso)--> Running ruff format $(ccend)"
	$(PYTHON_ENV) -m ruff format steps/calib_dedup/ steps/fastp/ steps/vidjil/ steps/igblast/ steps/cdr3nt_error_corrector/

integration-tests: venv ## >> run tests for all workflows via pytest and pytest-workflow tool
	@echo ""
	@echo "$(ccso)--> Running workflow tests $(ccend)"
	$(PYTHON_ENV) -m pytest tests/ -vv

unit-tests: venv ## >> run tests for all steps via pytest tool
	@echo ""
	@echo "$(ccso)--> Running steps tests $(ccend)"
#	$(PYTHON_ENV) -m pytest steps/calib_dedup/unit_tests -vv
#	$(PYTHON_ENV) -m pytest steps/fastp/unit_tests -vv
#	$(PYTHON_ENV) -m pytest steps/vidjil/unit_tests -vv
	$(PYTHON_ENV) -m pytest steps/igblast/unit_tests -vv
	$(PYTHON_ENV) -m pytest steps/cdr3nt_error_corrector/unit_tests -vv

tests: venv ##@main >> run integration and unit tests
	@echo ""
	@echo "$(ccso)--> Running integration and unit tests $(ccend)"
	$(MAKE) unit-tests
	$(MAKE) integration-tests

clean: ## >> remove docker images, python environment and nextflow build files
	@echo ""
	@echo "$(ccso)--> Removing temporary files and images $(ccend)"
	docker rmi -f downloader \
		calib-dedup-tool calib-dedup-image \
		fastp-tool fastp-image \
		vidjil-tool vidjil-image \
		igblast-tool igblast-image \
		cdr3nt-error-corrector-tool cdr3nt-error-corrector-image
	rm -rf $(VIRTUAL_ENV) \
		.nextflow.log* work .nextflow nextflow \
		/tmp/pytest_workflow_*

build-ref-images:
	@echo ""
	@echo "$(ccso)--> Build images of reference generators $(ccend)"
	docker build --target build-ref -t igblast-$(BUILD_REF_STAGE) steps/igblast/

build-igblast-ref-major: ## >> build an archive with igblast vdj reference with only major allele (*01)
	@echo ""
	@echo "$(ccso)--> Build a vdj reference with all alleles (*01) for igblast $(ccend)"
	docker run --rm -v ./steps/igblast:/work igblast-$(BUILD_REF_STAGE) -o /work/igblast.reference.major_allele.tar.gz

build-igblast-ref-all: ## >> build an archive with igblast vdj reference with all alleles
	@echo ""
	@echo "$(ccso)--> Build a vdj reference with all alleles (*01, *02, etc.) for igblast $(ccend)"
	docker run --rm -v ./steps/igblast:/work igblast-$(BUILD_REF_STAGE) -a -o /work/igblast.reference.all_alleles.tar.gz

build-vidjil-ref: ## >> build an archive with vidjil reference
	@echo ""
	@echo "$(ccso)--> Build a vdj reference for vidjil $(ccend)"
	bash steps/vidjil/build_ref.sh
	mv /tmp/vidjil.germline.tar.gz steps/vidjil/

build-olga-models: ## >> build an archive with olga models
	@echo ""
	@echo "$(ccso)--> Build olga models for cdr3nt-error-corrector $(ccend)"
	bash steps/cdr3nt_error_corrector/build_ref.sh
	mv /tmp/olga-models.tar.gz steps/cdr3nt_error_corrector/

build-ref: ##@main >> build all references
	@echo ""
	@echo "$(ccso)--> Build all references $(ccend)"
	$(MAKE) build-ref-images
	$(MAKE) build-igblast-ref-major
	$(MAKE) build-igblast-ref-all
	$(MAKE) build-vidjil-ref
	$(MAKE) build-olga-models

build: ##@main >> build docker images, the virtual environment and install requirements
	@echo ""
	@echo "$(ccso)--> Build $(ccend)"
	$(MAKE) clean
	$(MAKE) install-nextflow
	$(MAKE) build-images STAGE=image
	$(MAKE) build-images STAGE=tool
	$(MAKE) update

venv: $(VIRTUAL_ENV)

$(VIRTUAL_ENV): ## >> setup the virtual environment
	@echo "$(ccso)--> Install and setup virtualenv $(ccend)"
	$(PYTHON_SYS) -m pip install --upgrade pip
	$(PYTHON_SYS) -m pip install virtualenv
	virtualenv $(VIRTUAL_ENV)

update: venv ## >> update requirements.txt inside the virtual environment
	@echo "$(ccso)--> Updating packages $(ccend)"
	$(PYTHON_ENV) -m pip install -r ./steps/calib_dedup/requirements.txt
	$(PYTHON_ENV) -m pip install -r ./steps/igblast/requirements.txt
	$(PYTHON_ENV) -m pip install -r ./steps/cdr3nt_error_corrector/requirements.txt
	$(PYTHON_ENV) -m pip install pytest==8.1.1 pytest-workflow==2.1.0 ruff==0.4.2 mypy==1.10.0

install-python: ## >> install a python
	curl -fL http://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz > /tmp/python.tgz
	mkdir -p /tmp/python
	tar xzvf /tmp/python.tgz --one-top-level=/tmp/python --strip-component 1
	cd /tmp/python && \
		./configure --enable-optimizations && \
		make -j 8 && sudo make altinstall
	sudo apt install ${PYTHON_SYS}-distutils
	sudo rm -rf /tmp/python.tgz /tmp/python

install-docker: ## >> install a docker
	curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
	sudo add-apt-repository "deb [arch=${ARCHITECTURE}] https://download.docker.com/linux/ubuntu $(shell lsb_release -cs) stable"
	sudo apt-get update
	sudo apt-get -y -o Dpkg::Options::="--force-confnew" install docker-ce

install-java: ## >> install a java
	curl -fL https://download.oracle.com/java/${JAVA_VERSION}/latest/jdk-${JAVA_VERSION}_linux-x64_bin.tar.gz > /tmp/java.tar.gz
	mkdir -p /tmp/java-${JAVA_VERSION}
	tar xzvf /tmp/java.tar.gz --one-top-level=/tmp/java-${JAVA_VERSION} --strip-component 1
	sudo mv /tmp/java-${JAVA_VERSION} /usr/local/bin/
	@echo "export JAVA_HOME=/usr/local/bin/java-${JAVA_VERSION}" >> ~/.bashrc
	@echo 'export PATH=$$JAVA_HOME/bin:$$PATH' >> ~/.bashrc
	@java -version
	rm /tmp/java.tar.gz

install-nextflow: ## >> install a NextFlow
	@java --version
	curl -s https://get.nextflow.io | bash

build-images:
	@docker version
	docker build -t downloader steps/downloader
	docker build --target $(STAGE) -t calib-dedup-$(STAGE) steps/calib_dedup
	docker build --target $(STAGE) -t fastp-$(STAGE) steps/fastp
	docker build --target $(STAGE) -t vidjil-$(STAGE) steps/vidjil
	docker build --target $(STAGE) -t igblast-$(STAGE) steps/igblast
	docker build --target $(STAGE) -t cdr3nt-error-corrector-$(STAGE) steps/cdr3nt_error_corrector

install: ## Install and check dependencies
	$(MAKE) install-nextflow
	$(MAKE) build-ref
	$(MAKE) build-images STAGE=image

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
