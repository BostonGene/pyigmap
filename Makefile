install-java: ## >> Install Java
	curl -fL https://download.oracle.com/java/${JAVA_VERSION}/latest/jdk-${JAVA_VERSION}_linux-x64_bin.tar.gz -o /tmp/java.tar.gz
	mkdir -p /usr/local/bin/java-${JAVA_VERSION}
	tar --strip-components=1 -xzvf /tmp/java.tar.gz -C /usr/local/bin/java-${JAVA_VERSION}
	rm /tmp/java.tar.gz
	echo "export JAVA_HOME=/usr/local/bin/java-${JAVA_VERSION}" | tee -a ~/.bashrc ~/.zshrc
	echo 'export PATH=$$JAVA_HOME/bin:$$PATH' | tee -a ~/.bashrc ~/.zshrc
	export JAVA_HOME=/usr/local/bin/java-${JAVA_VERSION}
	export PATH=$$JAVA_HOME/bin:$$PATH