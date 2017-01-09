##
## JFX applications
##

define fun_jfx1
	mkdir -p webstart/tmp/com/github/lindenb/jvarkit/tools/jfx/$(dir $(1))
	xsltproc \
		--xinclude \
		-o webstart/tmp/com/github/lindenb/jvarkit/tools/jfx/$(1).fxml \
		src/main/resources/xsl/webstart2jfxml.xsl \
		src/main/java/com/github/lindenb/jvarkit/tools/jfx/$(1).xml
	xsltproc \
		--xinclude \
		--stringparam name "$(notdir $(1))" \
		--stringparam codebase "${webstart.base}" \
		--stringparam mainclass "com.github.lindenb.jvarkit.tools.jfx.$(subst /,.,$(1))" \
		-o webstart/$(notdir $(1)).jnlp \
		src/main/resources/xsl/webstart2jnlp.xsl \
		src/main/java/com/github/lindenb/jvarkit/tools/jfx/$(1).xml
	xsltproc \
		--xinclude \
		--stringparam name "$(notdir $(1))" \
		src/main/resources/xsl/webstart2html.xsl \
		src/main/java/com/github/lindenb/jvarkit/tools/jfx/$(1).xml >> webstart/index.html ;
endef

define sign_jfx1
	jarsigner -tsa http://timestamp.digicert.com \
		-keystore .secret.keystore \
		-keypass "$(if ${keytool.keypass},${keytool.keypass},KEYTOOLPASS)" \
		-storepass "$(if ${keytool.storepass},${keytool.storepass},KEYTOOLSTOREPASS)" \
		"$(1)" secret ;
endef

ifneq (${gatk.jar},)	

PICARDJFX=$(addprefix picardjfx/,FilterVcfJfx GatherVcfsJfx FindMendelianViolationsJfx)
GATKJFX=$(addprefix gatkjfx/,SelectVariantsJfx CombineVariantsJfx)

test-webstart: compile-webstart 
	java -cp webstart/gatkjfx.jar com.github.lindenb.jvarkit.tools.jfx.gatkjfx.CombineVariantsJfx

scp-webstart: compile-webstart
	scp -r webstart/* "${webstart.remotedir}"

compile-webstart : .secret.keystore ${lib.dir}/picard.jar \
	$(addprefix src/main/java/com/github/lindenb/jvarkit/tools/jfx/, $(addsuffix .java,${PICARDJFX} ${GATKJFX}) $(addsuffix .xml,${PICARDJFX} ${GATKJFX}))
	mkdir -p webstart/tmp
	echo "<html><body><table><tr><th>Name</th><th>Description</th></tr>" > webstart/index.html
	# compile GATK tools
	mkdir -p $(dir $@)
	unzip -o "${gatk.jar}"  -d webstart/tmp
	$(foreach P,${GATKJFX},$(call fun_jfx1,$P))
	javac -d webstart/tmp -cp "${gatk.jar}" -sourcepath  src/main/java $(addprefix src/main/java/com/github/lindenb/jvarkit/tools/jfx/, $(addsuffix .java,${GATKJFX}))
	cp ./src/main/java/com/github/lindenb/jvarkit/jfx/components/FilesChooserPane.fxml \
		webstart/tmp/com/github/lindenb/jvarkit/jfx/components/
	cp ./src/main/java/com/github/lindenb/jvarkit/jfx/components/FileChooserPane.fxml \
		webstart/tmp/com/github/lindenb/jvarkit/jfx/components/	
	jar cvf webstart/gatkjfx.jar -C webstart/tmp .
	rm -rf webstart/tmp
	# compile picard tools
	mkdir -p webstart/tmp
	$(foreach P,${PICARDJFX},$(call fun_jfx1,$P))
	javac -d webstart/tmp -cp ${lib.dir}/picard.jar -sourcepath  src/main/java $(addprefix src/main/java/com/github/lindenb/jvarkit/tools/jfx/, $(addsuffix .java,${PICARDJFX}))
	cp ./src/main/java/com/github/lindenb/jvarkit/jfx/components/FilesChooserPane.fxml \
		webstart/tmp/com/github/lindenb/jvarkit/jfx/components/
	cp ./src/main/java/com/github/lindenb/jvarkit/jfx/components/FileChooserPane.fxml \
		webstart/tmp/com/github/lindenb/jvarkit/jfx/components/
	jar cvf webstart/picardjfx.jar -C webstart/tmp .
	rm -rf webstart/tmp
	# sign jars
	cp ${lib.dir}/picard.jar webstart/
	$(call sign_jfx1,webstart/picardjfx.jar)
	$(call sign_jfx1,webstart/gatkjfx.jar)
	$(call sign_jfx1,webstart/picard.jar)
	echo "</table></body></html>" >> webstart/index.html
	chmod 755 webstart/*.html webstart/*.jar webstart/*.jnlp 
	

${lib.dir}/picard.jar:
	mkdir -p $(dir $@)
	wget -O "$@" "https://github.com/broadinstitute/picard/releases/download/2.8.1/$(notdir $@)"

endif






