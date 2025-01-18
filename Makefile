SHELL:=/bin/bash

JAVA_VERSION=21
ARCHITECTURE=amd64
BUILD_REF_STAGE=build-ref
ENGINE=docker
UV_BIN = '$(HOME)/.local/bin/uv'

ifeq ($(ENGINE),podman)
	PODMAN_PARAM=--format docker
	export USE_PODMAN=true
endif

# .ONESHELL:
DEFAULT_GOAL: install
.PHONY: help clean build build-ref tests integration-tests unit-tests mypy check format \
    install-python install-docker install-java install-podman

# Colors for echos 
ccend = $(shell tput sgr0)
ccbold = $(shell tput bold)
ccgreen = $(shell tput setaf 2)
ccso = $(shell tput smso)

mypy: venv ## >> run mypy type checker
	@echo ""
	@echo "$(ccso)--> Running mypy $(ccend)"
	$(UV_BIN) run mypy bin/calib_dedup/
	$(UV_BIN) run mypy bin/fastp/
	$(UV_BIN) run mypy bin/vidjil/
	$(UV_BIN) run mypy bin/igblast/
	$(UV_BIN) run mypy bin/cdr3nt_error_corrector/

check: ## >> run ruff linter
	@echo ""
	@echo "$(ccso)--> Running ruff check $(ccend)"
	$(UV_BIN) run ruff check bin/calib_dedup/ bin/fastp/ bin/vidjil/ bin/igblast/ bin/cdr3nt_error_corrector/

format: ## >> run ruff formatter
	@echo ""
	@echo "$(ccso)--> Running ruff format $(ccend)"
	$(UV_BIN) run ruff format bin/calib_dedup/ bin/fastp/ bin/vidjil/ bin/igblast/ bin/cdr3nt_error_corrector/

integration-tests: ## >> run tests for all workflows via pytest and pytest-workflow tool
	@echo ""
	@echo "$(ccso)--> Running workflow tests $(ccend)"
	$(UV_BIN) run pytest tests -vv --kwdof --tag integration-tests

unit-tests: ## >> run tests for all steps via pytest tool
	@echo ""
	@echo "$(ccso)--> Running steps tests $(ccend)"
	$(UV_BIN) run pytest bin/pyumi/unit_tests -vv
#	$(UV_BIN) run pytest bin/calib_dedup/unit_tests -vv
	$(UV_BIN) run pytest tests -vv --kwdof --tag unit-tests
	$(UV_BIN) run pytest bin/cdr3nt_error_corrector/unit_tests -vv

tests: ##@main >> run integration and unit tests
	@echo ""
	@echo "$(ccso)--> Running integration and unit tests $(ccend)"
	$(MAKE) unit-tests
	$(MAKE) integration-tests

clean: ## >> remove docker images, python environment and nextflow build files
	@echo ""
	@echo "$(ccso)--> Removing temporary files and images $(ccend)"
	$(ENGINE) rmi -f downloader-image \
		pyumi-tool pyumi-image \
		calib_dedup-tool calib_dedup-image \
		fastp-tool fastp-image \
		vidjil-tool vidjil-image \
		igblast-tool igblast-image \
		cdr3nt_error_corrector-tool cdr3nt_error_corrector-image
	rm -rf $(VIRTUAL_ENV) \
		.nextflow.log* work .nextflow nextflow

build-ref-image:
	@echo ""
	@echo "$(ccso)--> Build images of reference generators $(ccend)"
	$(ENGINE) build --target build-ref -t $(STEP)-$(BUILD_REF_STAGE) $(PODMAN_PARAM) bin/$(STEP)/

build-igblast-ref-major: ## >> build an archive with igblast vdj reference with only major allele (*01)
	@echo ""
	@echo "$(ccso)--> Build a vdj reference with all alleles (*01) for igblast $(ccend)"
	$(ENGINE) run --rm -v $(pwd)/bin/igblast:/work igblast-$(BUILD_REF_STAGE) -o /work/igblast.reference.major_allele.tar.gz

build-igblast-ref-all: ## >> build an archive with igblast vdj reference with all alleles
	@echo ""
	@echo "$(ccso)--> Build a vdj reference with all alleles (*01, *02, etc.) for igblast $(ccend)"
	$(ENGINE) run --rm -v $(pwd)/bin/igblast:/work igblast-$(BUILD_REF_STAGE) -a -o /work/igblast.reference.all_alleles.tar.gz

build-vidjil-ref: ## >> build an archive with vidjil reference
	@echo ""
	@echo "$(ccso)--> Build a vdj reference for vidjil $(ccend)"
	bash bin/vidjil/build_ref.sh
	mv /tmp/vidjil.germline.tar.gz bin/vidjil/

build-olga-models: ## >> build an archive with olga models
	@echo ""
	@echo "$(ccso)--> Build olga models for cdr3nt-error-corrector $(ccend)"
	bash bin/cdr3nt_error_corrector/build_ref.sh
	mv /tmp/olga-models.tar.gz bin/cdr3nt_error_corrector/

build-ref: ##@main >> build all references
	@echo ""
	@echo "$(ccso)--> Build all references $(ccend)"
	$(MAKE) build-ref-image STEP=igblast
	$(MAKE) build-igblast-ref-major
	$(MAKE) build-igblast-ref-all
	$(MAKE) build-vidjil-ref
	$(MAKE) build-olga-models

build: ##@main >> build docker images, the virtual environment and install requirements
	@echo ""
	@echo "$(ccso)--> Build $(ccend)"
	$(MAKE) clean
	$(MAKE) build-step-image STEP=downloader STAGE=image
	for step in pyumi calib_dedup fastp vidjil igblast cdr3nt_error_corrector ; do \
    	$(MAKE) build-step-image STEP=$$step STAGE=image ; \
    	$(MAKE) build-step-image STEP=$$step STAGE=tool ; \
	done
	$(MAKE) update
	# $(UV_BIN) run nf-core pipelines schema build --no-prompts
	chmod +x pyigmap

update: install-uv ## >> update requirements.txt inside the virtual environment
	@echo "$(ccso)--> Updating packages $(ccend)"
	$(UV_BIN) add -r bin/pyumi/requirements.txt
	$(UV_BIN) add -r bin/cdr3nt_error_corrector/requirements.txt
	$(UV_BIN) add "pytest>=8.1.1" "pytest-workflow>=2.1.0" "ruff>=0.4.2" "mypy>=1.10.0" "nf-core>=2.14.1" "nextflow>=24.04.2"

install-uv: ## >> Installs uv
	@echo ""
	@echo "Installing uv..."
	curl -LsSf https://astral.sh/uv/install.sh | sh

install-docker: ## >> install a docker
	curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
	sudo add-apt-repository "deb [arch=${ARCHITECTURE}] https://download.docker.com/linux/ubuntu $(shell lsb_release -cs) stable"
	sudo apt-get update
	sudo apt-get -y -o Dpkg::Options::="--force-confnew" install docker-ce

install-podman: ## >> install a Podman
	sudo apt-get update
	sudo apt-get -y install podman

install-java: ## >> Install Java
	curl -fL https://download.oracle.com/java/${JAVA_VERSION}/latest/jdk-${JAVA_VERSION}_linux-x64_bin.tar.gz -o /tmp/java.tar.gz
	sudo mkdir -p /usr/local/bin/java-${JAVA_VERSION}
	sudo tar --strip-components=1 -xzvf /tmp/java.tar.gz -C /usr/local/bin/java-${JAVA_VERSION}
	rm /tmp/java.tar.gz
	echo "export JAVA_HOME=/usr/local/bin/java-${JAVA_VERSION}" | tee -a ~/.bashrc ~/.zshrc
	echo 'export PATH=$$JAVA_HOME/bin:$$PATH' | tee -a ~/.bashrc ~/.zshrc
	export JAVA_HOME=/usr/local/bin/java-${JAVA_VERSION}
	export PATH=$$JAVA_HOME/bin:$$PATH

install-gh:
	curl -sS https://webi.sh/gh | sh

install-cliff: ## >> Install git-cliff tool
	@echo ""
	@echo "$(ccso)--> Installing git-cliff tool $(ccend)"
	curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
	$$HOME/.cargo/bin/cargo install --version 2.6.1 git-cliff

changelog: install-cliff ## >> Update CHANGELOG.md
	@echo ""
	@echo "$(ccso)--> Updating CHANGELOG.md $(ccend)"
	git-cliff --tag $(TAG) --bump -o CHANGELOG.md
	git add CHANGELOG.md && git commit -m "chore(release): prepare for $(TAG)"
	git show

release: install-gh changelog
	git tag $(TAG) && git push origin $(TAG)
	git push
	git cliff --strip all --latest | gh release create $(TAG) --notes-file -

build-step-image:
	@$(ENGINE) version
	$(ENGINE) build --target $(STAGE) -t $(STEP)-$(STAGE) $(PODMAN_PARAM) bin/$(STEP)

install: ## Install and check dependencies
	$(MAKE) update
	$(MAKE) build-ref
	$(MAKE) build-step-image STEP=downloader STAGE=image
	for step in pyumi calib_dedup fastp vidjil igblast cdr3nt_error_corrector ; do \
    	$(MAKE) build-step-image STEP=$$step STAGE=image ; \
	done
	chmod +x pyigmap

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