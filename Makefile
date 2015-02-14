#
# Starting moving from ANT to make
#
SHELL=/bin/bash
this.makefile=$(lastword $(MAKEFILE_LIST))
this.dir=$(dir $(realpath ${this.makefile}))

#need local settings ? create a file 'local.mk' in this directory
ifneq ($(realpath local.mk),)
include $(realpath local.mk)
endif

# proxy for curl, etc...
curl.proxy=$(if ${http.proxy.host}${http.proxy.port},-x "${http.proxy.host}:${http.proxy.port}",)
xjc.proxy=$(if ${http.proxy.host}${http.proxy.port}, -httpproxy "${http.proxy.host}:${http.proxy.port}" ,)

ANT?=ant
JAVAC?=${JAVA_HOME}/bin/javac
JAVA?=${JAVA_HOME}/bin/java
JAR?=${JAVA_HOME}/bin/jar
htsjdk.version?=1.128
htsjdk.home?=${this.dir}htsjdk-${htsjdk.version}
htsjdk.jars=$(addprefix ${htsjdk.home}/dist/,$(addsuffix .jar,commons-jexl-2.1.1 commons-logging-1.1.1 htsjdk-${htsjdk.version} snappy-java-1.0.3-rc3))
src.dir=${this.dir}src/main/java
generated.dir=${this.dir}src/main/generated-sources
tmp.dir=${this.dir}_tmp-${htsjdk.version}
tmp.mft=${tmp.dir}//META-INF/MANIFEST.MF
dist.dir?=${this.dir}dist-${htsjdk.version}
biostars.id=59647 86363 86480 84452 90204 94573 103303 106668 130456

## http://stackoverflow.com/questions/9551416
EMPTY :=
SPACE := $(EMPTY) $(EMPTY)

define compile-htsjdk-cmd

## 1 : target name
## 2 : qualified main class name
## 3 : other deps

$(1)  : ${htsjdk.jars} \
		${generated.dir}/java/com/github/lindenb/jvarkit/util/htsjdk/HtsjdkVersion.java \
		$(if $(3),$(3), $(addsuffix .java,$(addprefix ${src.dir}/,$(subst .,/,$(2)))) )
	mkdir -p ${tmp.dir}/META-INF ${dist.dir}
	cp src/main/resources/messages/messages.properties ${tmp.dir}
	${JAVAC} -d ${tmp.dir} -g -classpath "$$(subst $$(SPACE),:,$$(filter %.jar,$$^))" -sourcepath ${src.dir}:${generated.dir}/java $$(filter %.java,$$^)
	#create META-INF/MANIFEST.MF
	echo "Manifest-Version: 1.0" > ${tmp.mft}
	echo "Main-Class: $(2)" >> ${tmp.mft}
	echo "Class-Path: $$(filter %.jar,$$^) ${dist.dir}/$(1).jar" | fold -w 71 | awk '{printf("%s%s\n",(NR==1?"": " "),$$$$0);}' >>  ${tmp.mft}
	echo -n "Git-Hash: " >> ${tmp.mft}
	$$(if $$(realpath .git/refs/heads/master),cat $$(realpath .git/refs/heads/master), echo "undefined")  >> ${tmp.mft} 
	echo -n "Compile-Date: " >> ${tmp.mft}
	date +%Y-%m-%d:%H-%m-%S >> ${tmp.mft}
	#create jar
	${JAR} cfm ${dist.dir}/$(1).jar ${tmp.mft}  -C ${tmp.dir} .
	#create bash executable
	echo '#!/bin/bash' > ${dist.dir}/$(1)
	echo '${JAVA} -Xmx500m -cp "$$(subst $$(SPACE),:,$$(filter %.jar,$$^)):${dist.dir}/$(1).jar" $(2) $$*' > ${dist.dir}/$(1)
	chmod  ugo+rx ${dist.dir}/$(1)
	#cleanup
	rm -rf ${tmp.dir}

endef

#
# $1 :biostar post-id
# $2: other deps
#
define compile_biostar_cmd
$(call compile-htsjdk-cmd,biostar$(1),com.github.lindenb.jvarkit.tools.biostar.Biostar$(1),$(2))
endef

# 
# All executables
#
biostars: $(foreach B, ${biostars.id} , biostar$(B) )
APPS=vcfresetvcf sam2tsv vcffilterjs vcfgo biostars

.PHONY: all $(APPS) clean 

all: $(APPS)


$(eval $(call compile-htsjdk-cmd,sam2tsv,com.github.lindenb.jvarkit.tools.sam2tsv.Sam2Tsv,${src.dir}/com/github/lindenb/jvarkit/tools/sam2tsv/Sam2Tsv.java))
$(eval $(call compile-htsjdk-cmd,vcfresetvcf,com.github.lindenb.jvarkit.tools.misc.VcfRemoveGenotypeIfInVcf))
$(eval $(call compile-htsjdk-cmd,vcfgo,com.github.lindenb.jvarkit.tools.vcfgo.VcfGeneOntology))

$(eval $(foreach B, ${biostars.id} , $(call compile_biostar_cmd,$B)))

$(filter-out ${htsjdk.home}/dist/htsjdk-${htsjdk.version}.jar  ,${htsjdk.jars}) : ${htsjdk.home}/dist/htsjdk-${htsjdk.version}.jar 
	touch --no-create $@

${htsjdk.home}/dist/htsjdk-${htsjdk.version}.jar : ${htsjdk.home}/build.xml
	echo "Compiling htsjdk with $${JAVA_HOME} = ${JAVA_HOME}"
	(cd ${htsjdk.home} && ${ANT} )

${htsjdk.home}/build.xml : 
	mkdir -p $(dir ${htsjdk.home})
	rm -rf $(dir ${htsjdk.home})${htsjdk.version}.zip $(dir $@) 
	echo "Downloading HTSJDK ${htsjdk.version} with curl"
	curl  ${curl.proxy} -o $(dir ${htsjdk.home})${htsjdk.version}.zip -L "https://github.com/samtools/htsjdk/archive/${htsjdk.version}.zip"
	unzip $(dir ${htsjdk.home})${htsjdk.version}.zip -d $(dir ${htsjdk.home})
	find ${htsjdk.home} -exec touch '{}'  ';'
	rm -f $(dir ${htsjdk.home})${htsjdk.version}.zip

${generated.dir}/java/com/github/lindenb/jvarkit/util/htsjdk/HtsjdkVersion.java : ${htsjdk.home}/build.xml
	mkdir -p $(dir $@)
	echo "package com.github.lindenb.jvarkit.util.htsjdk;" > $@
	echo '@javax.annotation.Generated("jvarkit")' >> $@
	echo 'public class HtsjdkVersion{ private HtsjdkVersion(){}' >> $@
	echo 'public static String getVersion() {return "${htsjdk.version}";}' >> $@
	echo 'public static String getHome() {return "${htsjdk.home}";}' >> $@
	echo '}'  >> $@

## API EVS
src/main/generated-sources/java/edu/washington/gs/evs/package-info.java :
	mkdir -p ${generated.dir}/java
	${JAVA_HOME}/bin/xjc ${xjc.proxy} -d ${generated.dir}/java \
		-p edu.washington.gs.evs \
		"http://evs.gs.washington.edu/wsEVS/EVSDataQueryService?wsdl"

ifeq ($(realpath /home/lindenb/package/knime_2.11.1/plugins/org.knime.core_2.11.1.0045704/knime-core.jar),)

jvarkit.knime.version=1.0.0
## Knime plugin for jvarkit:
${dist.dir}/knime/eclipse.jvarkit-${jvarkit.knime.version}.jar : ${dist.dir}/knime/knime.htsjdk-${htsjdk.version}.jar 
	mkdir -p ${tmp.dir}/META-INF $(dir $@)
	echo "Manifest-Version: 1.0" > ${tmp.mft}
	echo "Bundle-ManifestVersion: 2" >> ${tmp.mft}
	echo "Bundle-Name:  jvarkit for knime" >> ${tmp.mft}
	echo "Bundle-SymbolicName:  com.github.lindenb.jvarkit" >> ${tmp.mft}
	echo "Bundle-Version:  ${jvarkit.knime.version}" >> ${tmp.mft}	
	echo "Bundle-ClassPath: jvarkit.jar" >>  ${tmp.mft}	
	echo "Bundle-Vendor: Pierre Lindenbaum" >>  ${tmp.mft}	
	echo "BRequire-Bundle: org.eclipse.core.runtime," >>  ${tmp.mft}	
	echo " org.knime.workbench.core," >>  ${tmp.mft}	
	echo " org.knime.workbench.repository," >>  ${tmp.mft}	
	echo " org.knime.base," >>  ${tmp.mft}	
	echo " htsjdk.knime;bundle-version=\"${htsjdk.version}\"" >>  ${tmp.mft}	
	echo "Bundle-ActivationPolicy: lazy" >>  ${tmp.mft}	
	echo "Export-Package: com.github.lindenb.jvarkit.knime4bio" >>  ${tmp.mft}	
	echo "Bundle-RequiredExecutionEnvironment: JavaSE-1.7" >>  ${tmp.mft}	
	echo "<?xml version='1.0'>" > ${tmp.dir}/plugin.xml
	

## Knime plugin for htsjdk
${dist.dir}/knime/knime.htsjdk-${htsjdk.version}.jar : ${htsjdk.jars}
	mkdir -p ${tmp.dir} $(dir ${tmp.mft}) $(dir $@)
	cp $(filter %.jar,$^) ${tmp.dir}
	echo "Manifest-Version: 1.0" > ${tmp.mft}
	echo "Bundle-ManifestVersion: 2" >> ${tmp.mft}
	echo "Bundle-Name: htsjdk  for KNIME" >> ${tmp.mft}
	echo "Bundle-SymbolicName: htsjdk.knime" >> ${tmp.mft}
	echo "Bundle-Version: ${htsjdk.version}" >> ${tmp.mft}
	echo -n "Bundle-ClassPath: " > ${tmp.dir}/package.list
	echo "$(notdir $(filter %.jar,$^))" | tr " " "," >> ${tmp.dir}/package.list
	cat ${tmp.dir}/package.list | fold -w 71 | awk '{printf("%s%s\n",(NR==1?"": " "),$$0);}' >> ${tmp.mft}
	echo -n "Export-Package: " > ${tmp.dir}/package.list
	$(foreach J,$(filter %.jar,$^), jar tf ${J} | grep -v '-' | grep '/$$' | sed 's%/$$%%' | tr "/" "." >> ${tmp.dir}/package.list ; )
	cat ${tmp.dir}/package.list | tr "\n" "," | sed 's/,$$//' | fold -w 71 | awk '{printf("%s%s\n",(NR==1?"": " "),$$0);}' >> ${tmp.mft}
	rm -f ${tmp.dir}/package.list
	echo "Bundle-RequiredExecutionEnvironment: JavaSE-1.7" >> ${tmp.mft}
	${JAR} cfvm $@ ${tmp.mft}  -C ${tmp.dir} .
	rm -rf ${tmp.dir}

endif

knime: ../xslt-sandbox/stylesheets/knime/knime2java.xsl src/knime/resources/model/plugin.xml
	xsltproc --stringparam base.dir src/knime/generated-sources/java  $^


clean:
	rm -rf ${dist.dir}


