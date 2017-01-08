##
## JFX applications
##

define fun_jfx1
	mkdir -p webstart/tmp/com/github/lindenb/jvarkit/tools/jfx/$(dir $(1))
	xsltproc \
		-o webstart/tmp/com/github/lindenb/jvarkit/tools/jfx/$(1).fxml \
		src/main/resources/xsl/webstart2jfxml.xsl \
		src/main/java/com/github/lindenb/jvarkit/tools/jfx/$(1).xml
	xsltproc \
		--stringparam name "$(notdir $(1))" \
		--stringparam codebase "." \
		--stringparam mainclass "com.github.lindenb.jvarkit.tools.jfx.$(subst /,.,$(1))" \
		-o webstart/$(notdir $(1)).jnlp \
		src/main/resources/xsl/webstart2jnlp.xsl \
		src/main/java/com/github/lindenb/jvarkit/tools/jfx/$(1).xml
	xsltproc \
		--stringparam name "$(notdir $(1))" \
		src/main/resources/xsl/webstart2jnlp.xsl >> webstart/index.html
endef

PICARDJFX=picardjfx/FilterVcfJfx

picardjfx : .secret.keystore ${lib.dir}/picard.jar $(addprefix src/main/java/com/github/lindenb/jvarkit/tools/jfx/, $(addsuffix .java,${PICARDJFX}) $(addsuffix .xml,${PICARDJFX}))
	mkdir -p webstart/tmp
	echo "<html><body><dl>" > webstart/index.html
	$(foreach P,${PICARDJFX},$(call fun_jfx1,$P))
	echo "</dl></body></html>" >> webstart/index.html
	javac -d webstart/tmp -cp ${lib.dir}/picard.jar -sourcepath  src/main/java $(filter %.java,$^)
	cp ./src/main/java/com/github/lindenb/jvarkit/jfx/components/FilesChooserPane.fxml \
		webstart/tmp/com/github/lindenb/jvarkit/jfx/components/
	cp ./src/main/java/com/github/lindenb/jvarkit/jfx/components/FileChooserPane.fxml \
		webstart/tmp/com/github/lindenb/jvarkit/jfx/components/
	jar cvf webstart/picardjfx.jar -C webstart/tmp .
	jarsigner  -keystore .secret.keystore \
		-keypass "$(if ${keytool.keypass},${keytool.keypass},KEYTOOLPASS)" \
		-storepass "$(if ${keytool.storepass},${keytool.storepass},KEYTOOLSTOREPASS)" \
		webstart/picardjfx.jar secret
	##rm -rf webstart/tmp
	

${lib.dir}/picard.jar:
	mkdir -p $(dir $@)
	wget -O "$@" "https://github.com/broadinstitute/picard/releases/download/2.8.1/$(notdir $@)"



