#
# Makefile for jvarkit
#
SHELL=/bin/bash
this.makefile=$(lastword $(MAKEFILE_LIST))
this.dir=$(dir $(realpath ${this.makefile}))


top:
	@echo "This  is the top target. Run 'make name-of-target' to build the desired target. Run 'make all' if you're Pierre Lindenbaum. See http://lindenb.github.io/jvarkit/  for a list of the available tools." 


#need local settings ? create a file 'local.mk' in this directory
ifneq ($(realpath local.mk),)
include $(realpath local.mk)
endif


# proxy for curl, etc...
curl.proxy=$(if ${http.proxy.host}${http.proxy.port},-x "${http.proxy.host}:${http.proxy.port}",)
xjc.proxy=$(if ${http.proxy.host}${http.proxy.port}, -httpproxy "${http.proxy.host}:${http.proxy.port}" ,)

ANT?=ant
JAVAC?=javac
JAVAH?=javah
JAVA?=java
JAVACC?=javacc
JAR?=jar


#
# include java libraries from Maven
#
include maven.mk


XJC?=$(JAVA) -cp "$(subst $(SPACE),:,${xjc_jars})" -DenableExternalEntityProcessing=true com.sun.tools.xjc.XJCFacade

src.dir=${this.dir}src/main/java
generated.dir=${this.dir}src/main/generated-sources
tmp.dir=${this.dir}_tmp-${htsjdk.version}
tmp.mft=${tmp.dir}/META-INF/MANIFEST.MF
export dist.dir?=${this.dir}dist

bigwig.version=20150429
bigwig.jar?=lib/BigWig.jar
bigwig.log4j.jar=$(dir ${bigwig.jar})/log4j-1.2.15.jar
bigwig.jars=${bigwig.jar} ${bigwig.log4j.jar}


jvarkit.package=com.github.lindenb.jvarkit


## http://stackoverflow.com/questions/9551416
EMPTY :=
SPACE := $(EMPTY) $(EMPTY)

define CRLF


endef

gatk.tools = org_broadinstitute_gatk_engine_CommandLineGATK \
		org_broadinstitute_gatk_tools_walkers_variantutils_CombineVariants \
		org_broadinstitute_gatk_tools_walkers_annotator_VariantAnnotator \
		org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration 


define make_gatk_code

${generated.dir}/resources/json/gatk/$$(lastword $$(subst _, ,$(1))).json :
	mkdir -p $$(dir $$@)
	curl -Lk $${curl.proxy} -o "$$(addsuffix .tmp,$$@)" "https://software.broadinstitute.org/gatk/gatkdocs/$(1).php.json"
	mv $$(addsuffix .tmp,$$@) $$@

${generated.dir}/java/com/github/lindenb/jvarkit/gatk/commands/Abstract$$(lastword $$(subst _, ,$(1))).java :${generated.dir}/resources/json/gatk/$$(lastword $$(subst _, ,$(1))).json gatkcodegen ${this.dir}src/main/resources/velocity/gatkcmd.vm
	mkdir -p $$(dir $$@) && ${JAVA} -jar ${dist.dir}/gatkcodegen.jar -T $$(word 3,$$^) $$< > $$@

endef




define compile-htsjdk-cmd

## 1 : target name
## 2 : qualified main class name
## 3 : other deps

$(1)  : ${htsjdk.jars} \
		$(addsuffix .java,$(addprefix ${src.dir}/,$(subst .,/,$(2)))) \
		$(3) 
	echo "### COMPILING $(1) ######"
	mkdir -p ${tmp.dir}/META-INF ${dist.dir}
	mkdir -p ${tmp.dir}/$(dir $(subst .,/,$(2)))
	cp -v "$(addsuffix .java,$(addprefix ${src.dir}/,$(subst .,/,$(2))))" "${tmp.dir}/$(dir $(subst .,/,$(2)))"
	echo '### Printing javac version : it should be OpenJdk 11. if Not, check your $$$${PATH}.'
	${JAVAC} -version 2>&1 > $$(addsuffix .jdkversion,$$@) && \
		cat $$(addsuffix .jdkversion,$$@) && \
		grep -E '11\.[0-9_]+\.[0-9_]+' $$(addsuffix .jdkversion,$$@) && \
		rm $$(addsuffix .jdkversion,$$@)
	#information about the java jre/sdk.
	echo -n "javac=" > "${tmp.dir}/META-INF/jdk.properties"
	which ${JAVAC} >> "${tmp.dir}/META-INF/jdk.properties"
	echo -n "jar=" >> "${tmp.dir}/META-INF/jdk.properties"
	which ${JAR} >> "${tmp.dir}/META-INF/jdk.properties"
	echo -n "java=" >> "${tmp.dir}/META-INF/jdk.properties"
	which ${JAVA} >> "${tmp.dir}/META-INF/jdk.properties"
	echo "classpath=$$(subst $$(SPACE),:,$$(filter %.jar,$$^))" >> "${tmp.dir}/META-INF/jdk.properties"
	echo -n "self=" >> "${tmp.dir}/META-INF/jdk.properties"
	echo "${dist.dir}/$(1)$(if ${standalone},-fat).jar" >> "${tmp.dir}/META-INF/jdk.properties"
	#compile
	${JAVAC} \
		-J-Djvarkit.libs.jars='$$(subst $$(SPACE),:,$$(filter %.jar,$$^))' \
		-J-Djvarkit.main.class='$(2)' \
		-J-Djvarkit.this.dir='${this.dir}' \
		-d ${tmp.dir} \
		-g -classpath "$$(subst $$(SPACE),:,$$(filter %.jar,$$^))" \
		-sourcepath ${src.dir}:${generated.dir}/java $$(filter %.java,$$^)
ifeq (${standalone},yes)
	$$(foreach J,$$(filter %.jar,$$^),unzip -o $${J} -d ${tmp.dir};) 
endif
	#create META-INF/MANIFEST.MF
	echo "Manifest-Version: 1.0" > ${tmp.mft}
	echo "Main-Class: $(2)" >> ${tmp.mft}
	echo "Htsjdk-Version: ${htsjdk.version}" >> ${tmp.mft}
ifneq (${standalone},yes)
	echo "Class-Path: $$(realpath $$(filter %.jar,$$^)) ${dist.dir}/$(1).jar" | fold -w 71 | awk '{printf("%s%s\n",(NR==1?"": " "),$$$$0);}' >>  ${tmp.mft}
endif
	echo -n "Git-Hash: " >> ${tmp.mft}
	$$(if $$(realpath .git/refs/heads/master),cat $$(realpath .git/refs/heads/master), echo "undefined")  >> ${tmp.mft} 
	echo -n "Compile-Date: " >> ${tmp.mft}
	date +%Y-%m-%d:%H-%m-%S >> ${tmp.mft}
	#create jar
	${JAR} cfm ${dist.dir}/$(1)$(if ${standalone},-fat).jar ${tmp.mft}  -C ${tmp.dir} .
	#create bash executable
	echo '#!/bin/bash' > ${dist.dir}/$(1)
	echo -n '${JAVA} -Djvarkit.log.name=$(1) -Dfile.encoding=UTF8 -Xmx500m $(if ${http.proxy.host},-Dhttp.proxyHost=${http.proxy.host})  $(if ${http.proxy.port},-Dhttp.proxyPort=${http.proxy.port}) ' >> ${dist.dir}/$(1)
ifeq (${standalone},yes)
	echo -n ' -jar "${dist.dir}/$(1)-fat.jar" '  >> ${dist.dir}/$(1)
else
	echo -n ' -cp "$$(subst $$(SPACE),:,$$(realpath $$(filter %.jar,$$^))):${dist.dir}/$(1).jar" $(2) '  >> ${dist.dir}/$(1)
endif
	echo '"$$$$@"' >> ${dist.dir}/$(1)
	chmod  ugo+rx ${dist.dir}/$(1)
	# generate markdown if needed
	-if [ "${TRAVIS}" != "true" ] && [ "${standalone}" != "yes" ] ; then \
		${JAVA} -Djvarkit.doc.dir=${this.dir}docs  -Djvarkit.src.dir=${this.dir}src/main/java -jar "${dist.dir}/$(1)$(if ${standalone},-fat).jar" --help --helpFormat make-doc  ;\
	fi
	#cleanup
	rm -rf "${tmp.dir}"


endef


##
## Creates a CGI script for java app
##
define compile-cgi-cmd

${dist.dir}/$(1).cgi : $(1)
	echo "#!/bin/bash" > $$@
	echo "PREFIX=$$(dirname $$0)" >> $$@
	echo "${JAVA} -Dfile.encoding=UTF8 -Xmx500m -Dprefs.file.xml=/var/www/cgi-bin/prefs.xml -jar $PREFIX/$(1).jar $$*" >> $$@

endef

#
# $1 :biostar post-id
# $2: other deps
#
define compile_biostar_cmd
$(call compile-htsjdk-cmd,biostar$(1),${jvarkit.package}.tools.biostar.Biostar$(1),$(2))
endef



gatk_apps:$(if ${gatk.jar},gatkwalkers,)
	

APPS= vcffixindels vcftail vcfhead vcfburdenfisherh vcfburdenfisherv vcfburdenmaf  \
	vcfmulti2oneallele vcfin vcffilterso vcffilterjs vcfburdensplitter vcfpredictions \
	vcfpolyx vcfburdenfiltergenes vcfbed vcfbedsetfilter gatk_apps vcftrio vcffamilies  groupbygene \
	 addlinearindextobed	allelefreqcalc	almostsortedvcf	backlocate	bam2fastq lowresbam2raster bam2raster	bam2svg \
	bam2xml bam2wig		bamcmpcoverage	bamindexreadnames	bamliftover	bamqueryreadnames \
	bamrenamechr	bamsnvwig	bamstats04	bamstats05 bamtreepack	batchigvpictures	bedliftover \
	bedrenamechr	biostar103303	biostar130456	biostar59647	biostar76892	biostar77288 \
	biostar77828	biostar78285	biostar78400	biostar81455	biostar84452	biostar84786	biostar86363 \
	biostar86480	biostar90204	msa2vcf	biostar95652 biostar139647	biostar145820 blast2sam reduceblast	blastmapannots \
	blastn2snp	buildwpontology	bwamemdigest	bwamemnop	cmpbams	cmpbamsandbuild	coveragenormalizer \
	downsamplevcf	evs2vcf	evs2xml	fastq2fasta kg2bed \
	fastqentropy	fastqgrep	fastqjs	fastqphred64to33	fastqrecordtreepack	fastqrevcomp	fastqshuffle \
	fastqsplitinterleaved	findallcoverageatposition	findavariation	findcorruptedfiles	findmyvirus	findnewsplicesites	fixvarscanmissingheader \
	fixvcf	fixvcfformat	fixvcfmissinggenotypes	gcanddepth	genomicjaspar	genscan	 \
	howmanybamdict	illuminadir	ilmnfastqstats	impactofduplicates \
	liftover2svg	mapuniprot	mergesplittedblast	ncbitaxonomy2xml metrics2xml ngsfilessummary	noemptyvcf \
	nozerovariationvcf	pademptyfastq	pubmeddump	pubmedorcidgraph pubmedfilterjs	referencetovcf	sam2json \
	sam2psl	sam2tsv	sam4weblogo	samclipindelfraction	samextractclip	samfindclippedregions	samfixcigar \
	samgrep	samjs	samshortinvert	samstats01	sigframe	sortvcfoninfo \
	sortvcfonref2	splitbam3	splitbytile	splitread \
	vcf2hilbert	vcf2ps	vcf2rdf	vcf2sql	vcf2xml	vcfannobam	 \
	vcfbiomart	vcfcadd	vcfcmppred	vcfcomm	vcfcompare	vcfcomparegt \
	vcfconcat	vcfcutsamples	 \
	vcfjaspar	vcfliftover	vcfmerge	vcfmulti2one \
	vcfrebase	vcfregistry.cgi	vcfregulomedb	vcfrenamechr	vcfrenamesamples \
	vcfresetvcf	vcfsetdict	vcfmakedict vcfshuffle	vcfsimulator	vcfstats vcfcombinetwosnvs vcfstripannot \
	vcftabixml	vcftreepack	 vcfvcf	worldmapgenome \
	uniprotfilterjs skipxmlelements vcfensemblvep vcfgroupbypop bamtile xcontaminations \
	biostar3654  bioalcidae bioalcidaejdk vcfburden  vcfreplacetag vcfindextabix \
	vcfpeekvcf vcfgetvariantbyindex  vcfmulti2oneinfo bedindextabix vcf2bam vcffilterxpath \
	biostar140111 pcrclipreads  extendrefwithreads pcrslicereads samjmx vcfjmx gtf2xml sortsamrefname biostar154220 \
	biostar160470 biostar165777 blastfilterjs vcfcomparecallers bamclip2insertion localrealignreads biostar170742 biostar172515 \
	biostar173114 samslop biostar175929 vcfcalledwithanothermethod biostar178713 \
	vcfremovegenotypejs vcfgenesplitter bamstats02 bamstats02view sammaskalignedbases biostar105754 gff2kg \
	bam2sql vcfinjectpedigree vcfburdenrscriptv vcffilternotinpedigree vcfderby01 vcf2zip pubmedgender pubmedmap vcfdoest splitvcf \
	forkvcf gbrowserhtml bim2vcf queue2make concatsam samreadlengthdistribution biostar214299 \
	vcfmovefilterstoinfo gatkcodegen cmpbams4 vcfeigen01 biostar234081 biostar234230  vcfgnomad vcf2svg mergeblastxml \
	vcfannotwithbeacon commbams samscansplitreads samretrieveseqandqual pubmedcodinglang  biostar251649 samcolortag vcf2table \
	variantsinwindow  knime2txt lumpyvcf2circos vcfucsc xsltstream vcfloopovergenes vcffilterjdk samjdk vcfnocall2homref \
	vcfserver tviewserver vcftrap prettysam vcfremoveunusedalt lumpysort samaddpi goutils gb2gff \
	indexcov2vcf samcustomsortjdk


.PHONY: all tests $(APPS) clean download_all_maven library top    burden 





all: $(APPS)

burden: vcfburden vcfburdensplitter vcfburdenfisherh vcfburdenfisherv vcfburdenmaf  vcfburdenfiltergenes vcfinjectpedigree vcfburdenrscriptv vcffilternotinpedigree vcfderby01 vcfmovefilterstoinfo


${dist.dir}/testsng.jar: testsng
tests: ${testng.jars} ${dist.dir}/testsng.jar
	-${JAVA} \
		$(if ${http.proxy.host},-Dhttp.proxyHost=${http.proxy.host} -Dhttps.proxyHost=${http.proxy.host}) \
		$(if ${http.proxy.port},-Dhttp.proxyPort=${http.proxy.port} -Dhttps.proxyPort=${http.proxy.port}) \
		-cp "$(subst $(SPACE),:,$(filter %.jar,$^))" org.testng.TestNG -parallel false \
		-suitename "jvarkit" -testname "jvarkit" \
		-log 2 -d "test-output" -testjar ${dist.dir}/testsng.jar
	rm -vf ${dist.dir}/testsng.jar

tests2: ${testng.jars} ${htsjdk.jars}  ${httpclient.libs} api.ncbi.gb  ${bigwig.jars}  ${mysql.jar} ${jetty.jars} ${common.math3.libs}  ${gson.jar} ${berkeleydb.jar}  ${jaxb.jars}
	rm -rf "${tmp.dir}"
	mkdir -p "${tmp.dir}"
	${JAVAC} -d ${tmp.dir} -cp "$(subst $(SPACE),:,$(filter %.jar,$^))" -sourcepath ${generated.dir}/java:src/test/java:src/main/java `find src/test/java -type f -name "*.java" | grep -i -v jfx | grep -v SimplePlotTest`
	${JAVA} \
		-Dhttp.nonProxyHosts="localhost" \
		$(if ${http.proxy.host},-Dhttp.proxyHost=${http.proxy.host} -Dhttps.proxyHost=${http.proxy.host}) \
		$(if ${http.proxy.port},-Dhttp.proxyPort=${http.proxy.port} -Dhttps.proxyPort=${http.proxy.port}) \
		-cp "$(subst $(SPACE),:,$(filter %.jar,$^)):${tmp.dir}"  org.testng.TestNG -parallel false -d "test-output" ./src/test/resources/testng.xml
	rm -rf "${tmp.dir}"


#bigwig
$(eval $(call compile-htsjdk-cmd,vcfbigwig,${jvarkit.package}.tools.vcfbigwig.VCFBigWig,${jcommander.jar} ${bigwig.jars} ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,vcfensemblreg,	${jvarkit.package}.tools.ensemblreg.VcfEnsemblReg,${bigwig.jars} ${jcommander.jar}))
$(eval $(call compile_biostar_cmd,105754,${bigwig.jar} ${jcommander.jar} ))
# common math
$(eval $(call compile-htsjdk-cmd,naivecnvdetector,${jvarkit.package}.tools.structvar.NaiveCnvDetector,${jcommander.jar} ${common.math3.libs} ${bigwig.jars}))
#berkeley
$(eval $(call compile-htsjdk-cmd,vcfphylotree,${jvarkit.package}.tools.phylo.VcfPhyloTree,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,ngsfilesscanner,${jvarkit.package}.tools.ngsfiles.NgsFilesScanner,${jcommander.jar} ${berkeleydb.jar}))
$(eval $(call compile_biostar_cmd,92368,${jcommander.jar} ${berkeleydb.jar}))
#mysql
$(eval $(call compile-htsjdk-cmd,vcfucsc,${jvarkit.package}.tools.vcfucsc.VcfUcsc,${jcommander.jar} ${mysql.jar}))

$(eval $(call compile-htsjdk-cmd,addlinearindextobed,${jvarkit.package}.tools.misc.AddLinearIndexToBed,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,allelefreqcalc,${jvarkit.package}.tools.misc.AlleleFrequencyCalculator,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bam2xml,${jvarkit.package}.tools.bam2xml.Bam2Xml,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,almostsortedvcf,${jvarkit.package}.tools.sortvcfonref.AlmostSortedVcf,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,backlocate,${jvarkit.package}.tools.backlocate.BackLocate,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bam2fastq,${jvarkit.package}.tools.fastq.BamToFastq,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bam2raster,${jvarkit.package}.tools.bam2graphics.Bam2Raster,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,lowresbam2raster,${jvarkit.package}.tools.bam2graphics.LowResBam2Raster,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bam2svg,${jvarkit.package}.tools.bam2svg.BamToSVG,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,sv2svg,${jvarkit.package}.tools.bam2svg.SvToSVG,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bam2wig,${jvarkit.package}.tools.bam2wig.Bam2Wig,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bamcmpcoverage,${jvarkit.package}.tools.misc.BamCmpCoverage,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,samaddpi,${jvarkit.package}.tools.misc.SamAddPI,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bamindexreadnames,${jvarkit.package}.tools.bamindexnames.BamIndexReadNames,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bamliftover,${jvarkit.package}.tools.liftover.BamLiftOver,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,bamqueryreadnames,${jvarkit.package}.tools.bamindexnames.BamQueryReadNames,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bamrenamechr,${jvarkit.package}.tools.misc.ConvertBamChromosomes,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bamstatsjfx,${jvarkit.package}.tools.bamstats01.BamStatsJfx,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,bamstats02,${jvarkit.package}.tools.bamstats01.BamStats02,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,bamstats02view,${jvarkit.package}.tools.bamstats01.BamStats02View,${jcommander.jar} bamstats02 ))
$(eval $(call compile-htsjdk-cmd,bamstats04,${jvarkit.package}.tools.bamstats04.BamStats04,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bamstats05,${jvarkit.package}.tools.bamstats04.BamStats05,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bamtreepack,${jvarkit.package}.tools.treepack.BamTreePack,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,batchigvpictures,${jvarkit.package}.tools.batchpicts.BatchIGVPictures,${jcommander.jar} copy.opendoc.odp.resources))
$(eval $(call compile-htsjdk-cmd,bedliftover,${jvarkit.package}.tools.liftover.BedLiftOver,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bedrenamechr,${jvarkit.package}.tools.misc.ConvertBedChromosomes,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,103303,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,130456,${jcommander.jar} ))
$(eval $(call compile_biostar_cmd,145820,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,59647,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,76892,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,77288,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,77828,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,78285,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,78400,${jcommander.jar} ${jaxb.jars}))
$(eval $(call compile_biostar_cmd,81455,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,84452,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,84786,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,86363,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,86480,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,90204,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,95652,${jcommander.jar}  ${jaxb.jars} api.ncbi.gb))
$(eval $(call compile_biostar_cmd,3654,${jcommander.jar} ${jaxb.jars} api.ncbi.insdseq api.ncbi.blast))
$(eval $(call compile_biostar_cmd,154220,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,140111,${jcommander.jar} ${jaxb.jars} api.ncbi.dbsnp.gt ${generated.dir}/java/gov/nih/nlm/ncbi/dbsnp/gt/package-info.java))
$(eval $(call compile_biostar_cmd,160470,${jcommander.jar} ${jaxb.jars} api.ncbi.blast ))
$(eval $(call compile_biostar_cmd,165777,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,170742,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,172515,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,173114,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,175929,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,178713,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,214299,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,234081,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,234230,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,251649,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,322664,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,332826,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,336589,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,352930,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,blast2sam,${jvarkit.package}.tools.blast2sam.BlastToSam,${jcommander.jar} ${jaxb.jars} api.ncbi.blast ))
$(eval $(call compile-htsjdk-cmd,gb2gff,${jvarkit.package}.tools.genbank.GenbankToGff3,${jcommander.jar} ${jaxb.jars} api.ncbi.gb ))
$(eval $(call compile-htsjdk-cmd,blastmapannots,${jvarkit.package}.tools.blastmapannots.BlastMapAnnotations,${jcommander.jar} ${jaxb.jars} api.ncbi.blast api.ncbi.gb ${generated.dir}/java/org/uniprot/package-info.java))
$(eval $(call compile-htsjdk-cmd,blastn2snp,${jvarkit.package}.tools.blast.BlastNToSnp,${jcommander.jar} ${jaxb.jars} api.ncbi.blast ))
$(eval $(call compile-htsjdk-cmd,reduceblast,${jvarkit.package}.tools.blast.ReduceBlast,${jcommander.jar} ${jaxb.jars} api.ncbi.blast ))
$(eval $(call compile-htsjdk-cmd,mergeblastxml,${jvarkit.package}.tools.blast.MergeBlastXml,${jcommander.jar} ${jaxb.jars} api.ncbi.blast))
$(eval $(call compile-htsjdk-cmd,buildwpontology,${jvarkit.package}.tools.misc.BuildWikipediaOntology,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bwamemdigest,${jvarkit.package}.tools.mem.BWAMemDigest,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bwamemnop,${jvarkit.package}.tools.mem.BWAMemNOp,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,cmpbams,${jvarkit.package}.tools.cmpbams.CompareBams,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,cmpbams4,${jvarkit.package}.tools.cmpbams.CompareBams4,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,commbams,${jvarkit.package}.tools.cmpbams.CommBams,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,cmpbamsandbuild,${jvarkit.package}.tools.cmpbams.CompareBamAndBuild,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,coveragenormalizer,${jvarkit.package}.tools.misc.CoverageNormalizer,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,downsamplevcf,${jvarkit.package}.tools.misc.DownSampleVcf,${jcommander.jar}))
###$(eval $(call compile-htsjdk-cmd,evs2bed,${jvarkit.package}.tools.evs2bed.DumpExomeVariantServerData,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,evs2vcf,${jvarkit.package}.tools.evs2bed.EvsToVcf,${jcommander.jar} ${jaxb.jars} api.evs))
$(eval $(call compile-htsjdk-cmd,evs2xml,${jvarkit.package}.tools.evs2bed.EvsDumpXml,${jcommander.jar} ${jaxb.jars} api.evs))
$(eval $(call compile-htsjdk-cmd,fastq2fasta,${jvarkit.package}.tools.misc.FastqToFasta,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,fastqentropy,${jvarkit.package}.tools.fastq.FastqEntropy,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,fastqgrep,${jvarkit.package}.tools.misc.FastqGrep,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,fastqjs,${jvarkit.package}.tools.fastq.FastqJavascript,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,fastqphred64to33,${jvarkit.package}.tools.fastq.ConvertPhred64toFastq33,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,fastqrecordtreepack,${jvarkit.package}.tools.treepack.FastqRecordTreePack,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,fastqrevcomp,${jvarkit.package}.tools.misc.FastqRevComp,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,fastqshuffle,${jvarkit.package}.tools.fastq.FastqShuffle,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,fastqsplitinterleaved,${jvarkit.package}.tools.fastq.FastqSplitInterleaved,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,findallcoverageatposition,${jvarkit.package}.tools.misc.FindAllCoverageAtPosition,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,findavariation,${jvarkit.package}.tools.misc.FindAVariation,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,findcorruptedfiles,${jvarkit.package}.tools.misc.FindCorruptedFiles,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,findmyvirus,${jvarkit.package}.tools.mem.FindMyVirus,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,findnewsplicesites,${jvarkit.package}.tools.rnaseq.FindNewSpliceSites,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,fixvarscanmissingheader,${jvarkit.package}.tools.misc.FixVarScanMissingVCFHeader,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,fixvcf,${jvarkit.package}.tools.misc.FixVCF,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,fixvcfformat,${jvarkit.package}.tools.misc.FixVcfFormat,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,fixvcfmissinggenotypes,${jvarkit.package}.tools.misc.FixVcfMissingGenotypes,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,gcanddepth,${jvarkit.package}.tools.misc.GcPercentAndDepth,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,concatsam,${jvarkit.package}.tools.misc.ConcatSam,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,genomicjaspar,${jvarkit.package}.tools.jaspar.GenomicJaspar,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,genscan,${jvarkit.package}.tools.genscan.GenScan,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,groupbygene,${jvarkit.package}.tools.groupbygene.GroupByGene,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,howmanybamdict,${jvarkit.package}.tools.misc.HowManyBamDict,${jcommander.jar}))
####$(eval $(call compile-htsjdk-cmd,idea20130924,${jvarkit.package}.tools.bwamempcr.Idea20130924))
$(eval $(call compile-htsjdk-cmd,illuminadir,${jvarkit.package}.tools.misc.IlluminaDirectory,${jcommander.jar}  ${gson.jar}))
$(eval $(call compile-htsjdk-cmd,ilmnfastqstats,${jvarkit.package}.tools.misc.IlluminaStatsFastq,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,impactofduplicates,${jvarkit.package}.tools.impactdup.ImpactOfDuplicates,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,kg2bed,${jvarkit.package}.tools.misc.KnownGenesToBed,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,liftover2svg,${jvarkit.package}.tools.liftover.LiftOverToSVG,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,mapuniprot,${jvarkit.package}.tools.misc.MapUniProtFeatures,${jcommander.jar} ${jaxb.jars} ${generated.dir}/java/org/uniprot/package-info.java))
$(eval $(call compile-htsjdk-cmd,mergesplittedblast,${jvarkit.package}.tools.blast.MergeSplittedBlast,${jcommander.jar} ${jaxb.jars} api.ncbi.blast))
$(eval $(call compile-htsjdk-cmd,metrics2xml,${jvarkit.package}.tools.metrics2xml.PicardMetricsToXML,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,ncbitaxonomy2xml,${jvarkit.package}.tools.misc.NcbiTaxonomyToXml,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,ngsfilessummary,${jvarkit.package}.tools.ngsfiles.NgsFilesSummary,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,noemptyvcf,${jvarkit.package}.tools.misc.NoEmptyVCF,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,nozerovariationvcf,${jvarkit.package}.tools.misc.NoZeroVariationVCF,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,pademptyfastq,${jvarkit.package}.tools.misc.PadEmptyFastq,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,pubmeddump,${jvarkit.package}.tools.pubmed.PubmedDump,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,pubmed404,${jvarkit.package}.tools.pubmed.Pubmed404,${jcommander.jar} ${httpclient.libs}))
$(eval $(call compile-htsjdk-cmd,ncbigenedump,${jvarkit.package}.tools.misc.NcbiGeneDump,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,pubmedcodinglang,${jvarkit.package}.tools.pubmed.PubmedCodingLanguages,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,pubmedgender,${jvarkit.package}.tools.pubmed.PubmedGender,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,pubmedmap,${jvarkit.package}.tools.pubmed.PubmedMap,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,pubmedgraph,${jvarkit.package}.tools.pubmed.PubmedGraph,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,pubmedorcidgraph,${jvarkit.package}.tools.pubmed.PubmedOrcidGraph,${jcommander.jar} ${berkeleydb.jar}))
$(eval $(call compile-htsjdk-cmd,pubmedauthorgraph,${jvarkit.package}.tools.pubmed.PubmedAuthorGraph,${jcommander.jar} ${berkeleydb.jar}))
$(eval $(call compile-htsjdk-cmd,pubmedfilterjs,${jvarkit.package}.tools.pubmed.PubmedFilterJS,${jcommander.jar} ${jaxb.jars} api.ncbi.pubmed))
$(eval $(call compile-htsjdk-cmd,referencetovcf,${jvarkit.package}.tools.misc.ReferenceToVCF,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,sam2json,${jvarkit.package}.tools.misc.SamToJson,${jcommander.jar} ${gson.jar} ))
$(eval $(call compile-htsjdk-cmd,sam2psl,${jvarkit.package}.tools.misc.SamToPsl,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,sam2tsv,${jvarkit.package}.tools.sam2tsv.Sam2Tsv,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,wescnvtview,${jvarkit.package}.tools.sam2tsv.WesCnvTView,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,prettysam,${jvarkit.package}.tools.sam2tsv.PrettySam,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,sam4weblogo,${jvarkit.package}.tools.sam4weblogo.SAM4WebLogo,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,samclipindelfraction,${jvarkit.package}.tools.misc.SamClipIndelFraction,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,samextractclip,${jvarkit.package}.tools.structvar.SamExtractClip,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,vcfstrechofgt,${jvarkit.package}.tools.structvar.VcfStretchOfGt,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,samscansplitreads,${jvarkit.package}.tools.structvar.SamScanSplitReads,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,samfindclippedregions,${jvarkit.package}.tools.structvar.SamFindClippedRegions,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,samfixcigar,${jvarkit.package}.tools.samfixcigar.SamFixCigar,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,samgrep,${jvarkit.package}.tools.samgrep.SamGrep,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,samjs,${jvarkit.package}.tools.samjs.SamJavascript,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,samcolortag,${jvarkit.package}.tools.samjs.SamColorTag,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,samshortinvert,${jvarkit.package}.tools.structvar.SamShortInvertion,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,samstats01,${jvarkit.package}.tools.bamstats01.BamStats01,${jcommander.jar} ))
##$(eval $(call compile-htsjdk-cmd,scanshortinvert,${jvarkit.package}.tools.mem.ScanShortInvert))
$(eval $(call compile-htsjdk-cmd,sigframe,${jvarkit.package}.tools.sigframe.SigFrame,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,sortvcfoninfo,${jvarkit.package}.tools.sortvcfonref.SortVcfOnInfo,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,sortvcfonref2,${jvarkit.package}.tools.sortvcfonref.SortVcfOnRef2,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,splitbam3,${jvarkit.package}.tools.splitbam.SplitBam3,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,splitbytile,${jvarkit.package}.tools.splitbytitle.SplitByTile,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,splitread,${jvarkit.package}.tools.splitread.SplitRead,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,tview,${jvarkit.package}.tools.tview.TViewCmd,${jcommander.jar}))
#$(eval $(call compile-cgi-cmd,tview.cgi))
$(eval $(call compile-htsjdk-cmd,vcf2hilbert,${jvarkit.package}.tools.hilbert.VcfToHilbert,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcf2ps,${jvarkit.package}.tools.misc.VcfToPostscript,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcf2svg,${jvarkit.package}.tools.misc.VcfToSvg,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcf2rdf,${jvarkit.package}.tools.vcf2rdf.VcfToRdf,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,vcf2sql,${jvarkit.package}.tools.vcf2sql.VcfToSql,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcf2xml,${jvarkit.package}.tools.vcf2xml.Vcf2Xml,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfannobam,${jvarkit.package}.tools.vcfannobam.VCFAnnoBam,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfbed,${jvarkit.package}.tools.vcfbed.VCFBed,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfbiomart,${jvarkit.package}.tools.ensembl.VcfBiomart,${jcommander.jar} ${httpclient.libs}))
$(eval $(call compile-htsjdk-cmd,vcfcadd,${jvarkit.package}.tools.misc.VcfCadd,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfallelebalance,${jvarkit.package}.tools.misc.VcfAlleleBalance,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfnocall2homref,${jvarkit.package}.tools.misc.VcfNoCallToHomRef,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfcmppred,${jvarkit.package}.tools.vcfcmp.VCFComparePredictions,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfcomm,${jvarkit.package}.tools.vcfcmp.VCFComm,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfcompare,${jvarkit.package}.tools.vcfcmp.VCFCompare,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfcomparegt,${jvarkit.package}.tools.vcfcmp.VCFCompareGT,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfconcat,${jvarkit.package}.tools.vcfconcat.VcfConcat,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcf2zip,${jvarkit.package}.tools.vcfconcat.VcfToZip,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,vcfcutsamples,${jvarkit.package}.tools.misc.VcfCutSamples,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfdas,${jvarkit.package}.tools.vcfdas.VcfDistributedAnnotationSystem,${jcommander.jar} ${jetty.jars} ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,vcffilterjs,${jvarkit.package}.tools.vcffilterjs.VCFFilterJS,${jcommander.jar}  ${gson.jar} ))
$(eval $(call compile-htsjdk-cmd,vcffilterso,${jvarkit.package}.tools.vcffilterso.VcfFilterSequenceOntology,${jcommander.jar} ${jaxb.jars} vcfpredictions))
$(eval $(call compile-htsjdk-cmd,vcffixindels,${jvarkit.package}.tools.vcffixindels.VCFFixIndels,${jcommander.jar}  ))
$(eval $(call compile-htsjdk-cmd,vcfgo,${jvarkit.package}.tools.vcfgo.VcfGeneOntology,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfhead,${jvarkit.package}.tools.misc.VcfHead,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcf2bed,${jvarkit.package}.tools.misc.VcfToBed,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,splitvcf,${jvarkit.package}.tools.misc.SplitVcf,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,forkvcf,${jvarkit.package}.tools.misc.ForkVcf,${jcommander.jar}  ))
$(eval $(call compile-htsjdk-cmd,vcfin,${jvarkit.package}.tools.vcfcmp.VcfIn,${jcommander.jar}  ))
$(eval $(call compile-htsjdk-cmd,vcfjaspar,${jvarkit.package}.tools.jaspar.VcfJaspar,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfliftover,${jvarkit.package}.tools.liftover.VcfLiftOver,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfmerge,${jvarkit.package}.tools.vcfmerge.VCFMerge,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfmulti2one,${jvarkit.package}.tools.onesamplevcf.VcfMultiToOne,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfpolyx,${jvarkit.package}.tools.misc.VCFPolyX,${jcommander.jar}   ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,vcfpredictions,${jvarkit.package}.tools.vcfannot.VCFPredictions,${jcommander.jar}  ))
$(eval $(call compile-htsjdk-cmd,vcfrebase,${jvarkit.package}.tools.vcfrebase.VcfRebase,${jcommander.jar}  ${jaxb.jars}) )
$(eval $(call compile-cgi-cmd,vcfregistry.cgi,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfregulomedb,${jvarkit.package}.tools.misc.VcfRegulomeDB,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfrenamechr,${jvarkit.package}.tools.misc.ConvertVcfChromosomes,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfrenamesamples,${jvarkit.package}.tools.misc.VcfRenameSamples,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfresetvcf,${jvarkit.package}.tools.misc.VcfRemoveGenotypeIfInVcf,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfremoveunusedalt,${jvarkit.package}.tools.misc.VcfRemoveUnusedAlt,${jcommander.jar}  ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,vcfremovegenotypejs,${jvarkit.package}.tools.misc.VcfRemoveGenotypeJs,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfsetdict,${jvarkit.package}.tools.misc.VcfSetSequenceDictionary,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfmakedict,${jvarkit.package}.tools.misc.VcfCreateDictionary,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfshuffle,${jvarkit.package}.tools.misc.VCFShuffle,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,vcfsimulator,${jvarkit.package}.tools.misc.VcfSimulator,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfstats,${jvarkit.package}.tools.vcfstats.VcfStats,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfstatsjfx,${jvarkit.package}.tools.vcfstats.VcfStatsJfx,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfcombinetwosnvs,${jvarkit.package}.tools.vcfannot.VCFCombineTwoSnvs,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfstripannot,${jvarkit.package}.tools.vcfstripannot.VCFStripAnnotations,${jcommander.jar} ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,vcftabixml,${jvarkit.package}.tools.vcftabixml.VCFTabixml,${jcommander.jar}  ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,vcftail,${jvarkit.package}.tools.misc.VcfTail,${jcommander.jar}  ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,vcftreepack,${jvarkit.package}.tools.treepack.VcfTreePack,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,vcftrio,${jvarkit.package}.tools.vcftrios.VCFTrios,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcffamilies,${jvarkit.package}.tools.vcftrios.VCFFamilies,${jcommander.jar} ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,vcfvcf,${jvarkit.package}.tools.vcfvcf.VcfVcf,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,worldmapgenome,${jvarkit.package}.tools.circular.WorldMapGenome,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,uniprotfilterjs,${jvarkit.package}.tools.misc.UniprotFilterJS,${jcommander.jar} ${jaxb.jars} ${generated.dir}/java/org/uniprot/package-info.java ))
$(eval $(call compile-htsjdk-cmd,blastfilterjs,${jvarkit.package}.tools.blast.BlastFilterJS,${jcommander.jar} ${jaxb.jars} api.ncbi.blast))
$(eval $(call compile-htsjdk-cmd,skipxmlelements,${jvarkit.package}.tools.misc.SkipXmlElements,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,minicaller,${jvarkit.package}.tools.calling.MiniCaller,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfcomparecallersonesample,${jvarkit.package}.tools.vcfcmp.VcfCompareCallersOneSample,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,samretrieveseqandqual,${jvarkit.package}.tools.misc.SamRetrieveSeqAndQual,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfensemblvep,${jvarkit.package}.tools.ensembl.VcfEnsemblVepRest,${jcommander.jar} ${httpclient.libs}))
$(eval $(call compile-htsjdk-cmd,vcfgroupbypop,${jvarkit.package}.tools.misc.VcfGroupByPopulation,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfcomparecallers,${jvarkit.package}.tools.vcfcmp.VcfCompareCallers,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bamtile,${jvarkit.package}.tools.misc.BamTile,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,xcontaminations,${jvarkit.package}.tools.xcontamination.XContaminations,${jcommander.jar}))
$(eval $(call compile_biostar_cmd,139647,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfburden,${jvarkit.package}.tools.misc.VcfBurden,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bioalcidae,${jvarkit.package}.tools.bioalcidae.BioAlcidae,${jcommander.jar} ${gson.jar} ${jaxb.jars} api.ncbi.blast api.ncbi.insdseq ${generated.dir}/java/gov/nih/nlm/ncbi/dbsnp/package-info.java  ))
$(eval $(call compile-htsjdk-cmd,bioalcidaejdk,${jvarkit.package}.tools.bioalcidae.BioAlcidaeJdk,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfbedsetfilter,${jvarkit.package}.tools.vcfbed.VCFBedSetFilter,${jcommander.jar}  ))
$(eval $(call compile-htsjdk-cmd,vcfreplacetag,${jvarkit.package}.tools.vcfstripannot.VCFReplaceTag,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfindextabix,${jvarkit.package}.tools.misc.VcfIndexTabix,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfpeekvcf,${jvarkit.package}.tools.vcfvcf.VcfPeekVcf,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfgetvariantbyindex,${jvarkit.package}.tools.misc.VcfGetVariantByIndex,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfmulti2oneallele,${jvarkit.package}.tools.misc.VcfMultiToOneAllele,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,vcfmulti2oneinfo,${jvarkit.package}.tools.misc.VcfMultiToOneInfo,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bedindextabix,${jvarkit.package}.tools.misc.BedIndexTabix,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcf2bam,${jvarkit.package}.tools.misc.VcfToBam,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcffilterxpath,${jvarkit.package}.tools.misc.VcfFilterXPath,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,pcrclipreads,${jvarkit.package}.tools.pcr.PcrClipReads,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bamslicebed,${jvarkit.package}.tools.pcr.BamSliceBed,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,extendrefwithreads,${jvarkit.package}.tools.extendref.ExtendReferenceWithReads,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,pcrslicereads,${jvarkit.package}.tools.pcr.PcrSliceReads,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,samjmx,${jvarkit.package}.tools.jmx.SamJmx,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfjmx,${jvarkit.package}.tools.jmx.VcfJmx,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,gtf2xml,${jvarkit.package}.tools.misc.Gtf2Xml,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,sortsamrefname,${jvarkit.package}.tools.misc.SortSamRefName,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bamclip2insertion,${jvarkit.package}.tools.misc.BamClipToInsertion,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,localrealignreads,${jvarkit.package}.tools.misc.LocalRealignReads,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,msa2vcf,${jvarkit.package}.tools.msa2vcf.MsaToVcf,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,samslop,${jvarkit.package}.tools.misc.SamSlop,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,projectserver,${jvarkit.package}.tools.server.ProjectServer,${jcommander.jar} ${jetty.jars}))
$(eval $(call compile-htsjdk-cmd,vcfcalledwithanothermethod,${jvarkit.package}.tools.misc.VcfCalledWithAnotherMethod,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfderby01,${jvarkit.package}.tools.burden.VcfDerby01,${jcommander.jar} ${derby.jars}))
$(eval $(call compile-htsjdk-cmd,vcfburdensplitter,${jvarkit.package}.tools.burden.VcfBurdenSplitter,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,casectrljfx,${jvarkit.package}.tools.burden.CaseControlJfx,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,indexcovjfx,${jvarkit.package}.tools.structvar.IndexCovJfx,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,indexcov2vcf,${jvarkit.package}.tools.structvar.IndexCovToVcf,${jcommander.jar} ${common.math3.libs}))
$(eval $(call compile-htsjdk-cmd,validatecnv,${jvarkit.package}.tools.structvar.ValidateCnv,${jcommander.jar} ${common.math3.libs}))
$(eval $(call compile-htsjdk-cmd,vcfburdensplitter2,${jvarkit.package}.tools.burden.VcfBurdenSplitter2,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfburdenfisherh,${jvarkit.package}.tools.burden.VcfBurdenFisherH,${jcommander.jar} ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,vcfburdenfisherv,${jvarkit.package}.tools.burden.VcfBurdenFisherV,${jcommander.jar} ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,vcfburdenrscriptv,${jvarkit.package}.tools.burden.VcfBurdenRscriptV,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfdoest,${jvarkit.package}.tools.burden.VcfDoest,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcffilternotinpedigree,${jvarkit.package}.tools.burden.VcfFilterNotInPedigree,${jcommander.jar} ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,vcfburdenmaf,${jvarkit.package}.tools.burden.VcfBurdenMAF,${jcommander.jar} ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,vcfburdenexac,${jvarkit.package}.tools.burden.VcfBurdenFilterExac,${jcommander.jar} ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,vcfgenesplitter,${jvarkit.package}.tools.misc.VcfGeneSplitter,${jcommander.jar} ))
#$(eval $(call compile-htsjdk-cmd,vcfsqltag,${jvarkit.package}.tools.sql.VcfSqlTag))
$(eval $(call compile-htsjdk-cmd,vcfburdenfiltergenes,${jvarkit.package}.tools.burden.VcfBurdenFilterGenes,${jcommander.jar}  ))
$(eval $(call compile-htsjdk-cmd,sammaskalignedbases,${jvarkit.package}.tools.misc.SamMaskAlignedBases,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,gff2kg,${jvarkit.package}.tools.misc.Gff2KnownGene,${jcommander.jar} ))
#$(eval $(call compile-htsjdk-cmd,miniassembly,${jvarkit.package}.tools.misc.MiniAssembly,))
$(eval $(call compile-htsjdk-cmd,haloplexparasite,${jvarkit.package}.tools.haloplex.HaloplexParasite,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,bam2sql,${jvarkit.package}.tools.misc.BamToSql,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfinjectpedigree,${jvarkit.package}.tools.burden.VcfInjectPedigree,${jcommander.jar} ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,gbrowserhtml,${jvarkit.package}.tools.misc.GBrowserHtml,${jcommander.jar}  ${gson.jar} copy.samtools.js))
$(eval $(call compile-htsjdk-cmd,bim2vcf,${jvarkit.package}.tools.misc.BimToVcf,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,samreadlengthdistribution,${jvarkit.package}.tools.misc.SamReadLengthDistribution,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,queue2make,${jvarkit.package}.tools.misc.QueueToMake,${jcommander.jar} ))
$(eval $(call compile-htsjdk-cmd,vcfmovefilterstoinfo,${jvarkit.package}.tools.burden.VcfMoveFiltersToInfo,${jcommander.jar} ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,gatkcodegen,${jvarkit.package}.tools.gatk.codegen.GATKCodeGenerator,${jcommander.jar} ${gson.jar} ${velocity.jars} ))
$(eval $(call compile-htsjdk-cmd,vcfeigen01,${jvarkit.package}.tools.vcfeigen.VcfEigen01,${jcommander.jar} ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,jfxngs,${jvarkit.package}.tools.vcfviewgui.JfxNgs,${jcommander.jar} jfxngs-resources))
$(eval $(call compile-htsjdk-cmd,ngsworkflow,${jvarkit.package}.tools.workflow.NgsWorkflow,${gson.jar} ${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfgnomad,${jvarkit.package}.tools.gnomad.VcfGnomad,${gson.jar} ${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfcomposite,${jvarkit.package}.tools.vcfcomposite.VCFComposite,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfannotwithbeacon,${jvarkit.package}.tools.ga4gh.VcfAnnotWithBeacon,${jcommander.jar} ${gson.jar} ${berkeleydb.jar} ${httpclient.libs} ))
$(eval $(call compile-htsjdk-cmd,vcf2table,${jvarkit.package}.tools.misc.VcfToTable,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,prettytable,${jvarkit.package}.tools.misc.PrettyTable,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,goutils,${jvarkit.package}.tools.misc.GoUtils,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,variantsinwindow,${jvarkit.package}.tools.misc.VariantsInWindow,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,casectrlcanvas,${jvarkit.package}.tools.burden.CaseControlCanvas,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,knime2txt,${jvarkit.package}.tools.misc.KnimeToText,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,lumpyvcf2circos,${jvarkit.package}.tools.lumpysv.LumpyVcfToCircos,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,lumpysort,${jvarkit.package}.tools.lumpysv.LumpySort,${jcommander.jar} ${berkeleydb.jar}))
$(eval $(call compile-htsjdk-cmd,fastgenotypegvcfs,${jvarkit.package}.tools.gvcf.FastGenotypeGVCFs,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,xsltstream,${jvarkit.package}.tools.misc.XsltStream,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfloopovergenes,${jvarkit.package}.tools.burden.VcfLoopOverGenes,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcffilterjdk,${jvarkit.package}.tools.vcffilterjs.VcfFilterJdk,${jcommander.jar} ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,samjdk,${jvarkit.package}.tools.samjs.SamJdk,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,samcustomsortjdk,${jvarkit.package}.tools.samjs.SamCustomSortJdk,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,optimizer,${jvarkit.package}.tools.optimizer.Optimizer,${jcommander.jar} ${gson.jar}))
##$(eval $(call compile-htsjdk-cmd,vcfamalgamation,${jvarkit.package}.tools.vcfamalgation.VcfXmlAmalgamation,${jcommander.jar} ${gson.jar}  ${bigwig.jars}))
$(eval $(call compile-htsjdk-cmd,vcfserver,${jvarkit.package}.tools.vcfserver.VcfServer,${jcommander.jar} ${jetty.jars}))
$(eval $(call compile-htsjdk-cmd,tviewserver,${jvarkit.package}.tools.tview.TViewServer,${jcommander.jar} ${jetty.jars}))
$(eval $(call compile-htsjdk-cmd,trapindexer,${jvarkit.package}.tools.trap.TrapIndexer,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcftrap,${jvarkit.package}.tools.trap.VcfTrap,${jcommander.jar}  ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,vcfepistatis01,${jvarkit.package}.tools.epistasis.VcfEpistatis01,${jcommander.jar}))
#$(eval $(call compile-htsjdk-cmd,subcloneit,${jvarkit.package}.tools.cloneit.SubCloneIt,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,igvreview,${jvarkit.package}.tools.igvreview.IgvReview,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,samtranslocations,${jvarkit.package}.tools.structvar.SamTranslocations,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfafinfofilter,${jvarkit.package}.tools.misc.VcfAfInfoFilter,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,faidxsplitter,${jvarkit.package}.tools.misc.FaidxSplitter,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfclusteredreadedge,${jvarkit.package}.tools.misc.VcfClusteredReadEdge,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfgapfrequent,${jvarkit.package}.tools.structvar.VcfGapFrequent,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,mergecnvnator,${jvarkit.package}.tools.structvar.MergeCnvNator,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,mergesv,${jvarkit.package}.tools.structvar.MergeStructuralVariants,${jcommander.jar}))


$(eval $(call compile-htsjdk-cmd,vcfburdengoenrichment,${jvarkit.package}.tools.burden.VcfBurdenGoEnrichment,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfoptimizeped4skat,${jvarkit.package}.tools.skat.VcfOptimizePedForSkat,${jcommander.jar} ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,vcfskatslidingwindow,${jvarkit.package}.tools.skat.VcfSkatSlidingWindow,${jcommander.jar}  ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,vcfskat,${jvarkit.package}.tools.skat.VcfSkat,${jcommander.jar}  ${jaxb.jars}))
##$(eval $(call compile-htsjdk-cmd,vcfspringfilter,${jvarkit.package}.tools.misc.VcfSpringFilter,${jcommander.jar} ${spring-beans.jars}))
##$(eval $(call compile-htsjdk-cmd,testsng,${jvarkit.package}.tools.tests.TestNg01,${testng.jars}  ${bigwig.jars} ${jaxb.jars}))
$(eval $(call compile-htsjdk-cmd,simpleplot,${jvarkit.package}.tools.misc.SimplePlot,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,cytoband2svg,${jvarkit.package}.tools.misc.CytobandToSvg,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,vcfancestralalleles,${jvarkit.package}.tools.onekgenomes.VcfAncestralAllele,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,bednonoverlappingset,${jvarkit.package}.tools.misc.BedNonOverlappingSet,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,wescnvsvg,${jvarkit.package}.tools.bam2svg.WesCnvSvg,${jcommander.jar}))
$(eval $(call compile-htsjdk-cmd,gexftr,${jvarkit.package}.tools.gephi.GexfTransformer,${jcommander.jar}))


$(lib.dir)/org/gephi/gephi.jar:
	test -f "${gephi_home}/bin/gephi"
	rm -rf "${tmp.dir}"
	mkdir -p "${tmp.dir}" "$(dir $@)"
	find "${gephi_home}/gephi/modules" \
		"${gephi_home}/platform/core" \
		"${gephi_home}/platform/lib" -type f -name "*.jar" | while read F; do echo "$$F" && unzip -o  "$$F" -d "${tmp.dir}" ; done
	jar cf $@ -C "${tmp.dir}" .
	rm -rf "${tmp.dir}"
	jar tvf $@ | head
	

$(eval $(call compile-htsjdk-cmd,gephicmd,${jvarkit.package}.tools.gephi.GephiCmd,${jcommander.jar} $(lib.dir)/org/gephi/gephi.jar ))

ij: ${dist.dir}/ij
${dist.dir}/ij : $(sort ${derby.jars} ${derby-tools.jar})
	mkdir -p $(dir $@)org.gephi.preview.api
	echo '#!/bin/bash' > @
	echo -n '${JAVA} -Dfile.encoding=UTF8  -cp "$(subst $(SPACE),:,$(filter %.jar,$^))" org.apache.derby.tools.ij $$*' >> $@
	chmod  ugo+rx $@

javacc: ${dist.dir}/javacc
${dist.dir}/javacc : ${javacc.jar}
	mkdir -p $(dir $@)
	echo '#!/bin/bash' > @
	echo -n '${JAVA} -Dfile.encoding=UTF8  -cp "$(subst $(SPACE),:,$(filter %.jar,$^))" javacc $$*' >> $@
	chmod  ugo+rx $@
	

GATKWALKERS_SRC=$(addsuffix .java,$(addprefix ${src.dir}/com/github/lindenb/jvarkit/tools/gatk/, variants/SoftClipAnnotator variants/GroupByVariants variants/GroupByGenotypes variants/EigenVariants variants/WindowVariants))

gatkwalkers:
	mkdir -p ${tmp.dir} ${dist.dir}
	${JAVAC} -d ${tmp.dir} -g -classpath ${gatk.jar} -sourcepath ${src.dir}:${generated.dir}/java ${GATKWALKERS_SRC}
	##$(foreach F,${GATKWALKERS_SRC}, cp ${F} ${tmp.dir}$(subst ${src.dir},,${F});)
	${JAR} cf ${dist.dir}/mygatk.jar -C ${tmp.dir} .
	rm -rf ${tmp.dir}
	echo '#!/bin/bash' > ${dist.dir}/mygatk
	echo -n '${JAVA} -Dfile.encoding=UTF8 -Xmx3g ' >> ${dist.dir}/mygatk
	echo -n ' -cp "$(realpath ${gatk.jar}):$(realpath ${dist.dir}/mygatk.jar)" org.broadinstitute.gatk.engine.CommandLineGATK '  >> ${dist.dir}/mygatk
	echo '$$*' >> ${dist.dir}/mygatk
	chmod  ugo+rx ${dist.dir}/mygatk


all-jnlp : $(addprefix ${dist.dir}/,$(addsuffix .jar, buildwpontology batchigvpictures)) ${htsjdk.jars} \
	 ./src/main/resources/jnlp/generic.jnlp .secret.keystore 
	mkdir -p ${tmp.dir}
	cp $(filter %.jar,$^) ${tmp.dir}
	$(foreach J, $(filter %.jar,$^) , jarsigner ${J} secret ; ) 
	sed -e 's/__MAIN_JAR__/buildwpontology/g' \
		-e 's/__TITLE__/BuildWikipediaOntology/g' \
		-e 's/__MAIN_CLASS__/${jvarkit.package}.tools.misc.BuildWikipediaOntology/g' \
		 ./src/main/resources/jnlp/generic.jnlp > ${tmp.dir}/buildwpontology.jnlp
	sed -e 's/__MAIN_JAR__/batchigvpictures/g' \
		-e 's/__TITLE__/BatchIGVPictures/g' \
		-e 's/__MAIN_CLASS__/${jvarkit.package}.tools.batchpicts.BatchIGVPictures/g' \
		 ./src/main/resources/jnlp/generic.jnlp > ${tmp.dir}/batchigvpictures.jnlp

.secret.keystore : 
	 -keytool -genkeypair -keystore $@ -alias secret \
	 		-keypass "$(if ${keytool.keypass},${keytool.keypass},KEYTOOLPASS)" \
	 		-storepass "$(if ${keytool.storepass},${keytool.storepass},KEYTOOLSTOREPASS)" \
	 		-dname "CN=Pierre Lindenbaum, OU=INSERM, O=INSERM, L=Nantes, ST=Nantes, C=Fr"

${generated.dir}/java/org/uniprot/package-info.java : api.uniprot
api.uniprot : ${xjc_jars}
	rm -rf ${generated.dir}/java/org/uniprot/
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java -p org.uniprot ${xjc.proxy} "https://www.uniprot.org/docs/uniprot.xsd" 


api.ncbi.pubmed : ${xjc_jars}
	rm -rf  ${generated.dir}/java/gov/nih/nlm/ncbi/pubmed/
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java  -p gov.nih.nlm.ncbi.pubmed -dtd ${xjc.proxy} "https://www.ncbi.nlm.nih.gov/corehtml/query/DTD/pubmed_100101.dtd"
	

api.evs: ${xjc_jars}
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java -p edu.washington.gs.evs ${xjc.proxy} -wsdl "http://evs.gs.washington.edu/wsEVS/EVSDataQueryService?wsdl"

api.ncbi.blast: ${xjc_jars}
	rm -rf  ${generated.dir}/gov/nih/nlm/ncbi/blast
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java  -p gov.nih.nlm.ncbi.blast -dtd ${xjc.proxy} "https://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd"

api.ncbi.esearch:
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java  -p gov.nih.nlm.ncbi.esearch -dtd ${xjc.proxy} "https://www.ncbi.nlm.nih.gov/corehtml/query/DTD/eSearch_020511.dtd"

api.ncbi.elink: ${xjc_jars}
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java  -p gov.nih.nlm.ncbi.elink -dtd ${xjc.proxy} "https://www.ncbi.nlm.nih.gov/corehtml/query/DTD/eLink_020511.dtd"


api.ncbi.gb: ${generated.dir}/java/gov/nih/nlm/ncbi/gb/ObjectFactory.java
${generated.dir}/java/gov/nih/nlm/ncbi/gb/ObjectFactory.java : ${xjc_jars}
	rm -rf  ${generated.dir}/gov/nih/nlm/ncbi/gb
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java  -p gov.nih.nlm.ncbi.gb -dtd ${xjc.proxy} "https://www.ncbi.nlm.nih.gov/dtd/NCBI_GBSeq.dtd"

api.ncbi.taxonomy: ${xjc_jars}
	rm -rf  ${generated.dir}/gov/nih/nlm/ncbi/taxonomy
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java  -p gov.nih.nlm.ncbi.taxonomy -dtd ${xjc.proxy} "https://www.ncbi.nlm.nih.gov/entrez/query/DTD/taxon.dtd"

api.ncbi.tseq: ${xjc_jars}
	rm -rf  ${generated.dir}/gov/nih/nlm/ncbi/tseq
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java  -p gov.nih.nlm.ncbi.tseq -dtd ${xjc.proxy} "https://www.ncbi.nlm.nih.gov/dtd/NCBI_TSeq.dtd"

api.ncbi.insdseq: ${xjc_jars}
	rm -rf  ${generated.dir}/gov/nih/nlm/ncbi/insdseq
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java  -p gov.nih.nlm.ncbi.insdseq -dtd ${xjc.proxy} "https://www.ncbi.nlm.nih.gov/dtd/INSD_INSDSeq.dtd"

api.ncbi.dbsnp.gt: ${this.dir}src/main/resources/xsd/ncbi/genoex_1_5.xsd ${jaxb.jars }
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java  -p gov.nih.nlm.ncbi.dbsnp.gt ${xjc.proxy} $<


${generated.dir}/java/gov/nih/nlm/ncbi/dbsnp/package-info.java : api.ncbi.dbsnp
	touch -c $@

api.ncbi.dbsnp: ${this.dir}src/main/resources/xsd/ncbi/docsum_3.4.xsd ${xjc_jars}
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java  -p gov.nih.nlm.ncbi.dbsnp ${xjc.proxy} $<


## API Ensembl
#api.ensembl.vep :
#	mkdir -p ${generated.dir}/java
#	${XJC} -d ${generated.dir}/java  -p org.ensembl.vep  ./src/main/resources/xsd/ensembl/vep.xsd


api.samfilter:
	mkdir -p ${src.dir}/java/com/github/lindenb/jvarkit/util/bio/samfilter && \
	${JAVACC} -OUTPUT_DIRECTORY=${src.dir}/com/github/lindenb/jvarkit/util/bio/samfilter \
		src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj

copy.opendoc.odp.resources : 
	mkdir -p ${tmp.dir}/META-INF/opendocument/odp
	cp src/main/resources/opendocument/thumbnail.png ${tmp.dir}/META-INF/opendocument/odp/
	cp $(addprefix src/main/resources/opendocument/odp/, meta.xml content.xml settings.xml styles.xml mimetype ) ${tmp.dir}/META-INF/opendocument/odp/

copy.samtools.js:
	mkdir -p ${tmp.dir}/META-INF/js
	$(foreach J,gbrowse.js hershey.js samtools.js com.github.lindenb.jvarkit.tools.misc.GBrowserHtml.js,cp src/main/js/${J} ${tmp.dir}/META-INF/js/ ; )

## GATK code generation using php.json URLs
gatk_code: $(foreach U,${gatk.tools},${generated.dir}/java/com/github/lindenb/jvarkit/gatk/commands/Abstract$(lastword $(subst _, ,$U)).java)
$(eval $(foreach U,${gatk.tools},$(call make_gatk_code,${U})))



## jvarkit-library (used in knime)
library: ${dist.dir}/jvarkit-${htsjdk.version}.jar
${dist.dir}/jvarkit-${htsjdk.version}.jar : ${htsjdk.jars} ${bigwig.jars} \
		${src.dir}/com/github/lindenb/jvarkit/util/Library.java
	mkdir -p ${tmp.dir}/META-INF $(dir $@)
	cp src/main/resources/messages/messages.properties ${tmp.dir}
	${JAVAC} -d ${tmp.dir} -g -classpath "$(subst $(SPACE),:,$(filter %.jar,$^))" -sourcepath ${src.dir}:${generated.dir}/java $(filter %.java,$^)
	${JAR} cf $@ -C ${tmp.dir} .
	rm -rf ${tmp.dir}

##
## Broad BigWig
## 
${bigwig.log4j.jar} : ${bigwig.jar}
	touch -c $@

${bigwig.jar} xx :
	echo "Downloading bigwig library for java."
	mkdir -p lib
	rm -rf bigwig.zip bigwig-${bigwig.version}
	curl -k ${curl.proxy} -o bigwig.zip -L "https://github.com/lindenb/bigwig/archive/${bigwig.version}.zip"
	unzip bigwig.zip
	rm -rf bigwig.zip 
	echo "Compiling, using apache ant"
	(cd  bigwig-${bigwig.version} && ${ANT} )
	mv bigwig-${bigwig.version}/lib/$(notdir ${bigwig.log4j.jar}) ${bigwig.log4j.jar}
	mv bigwig-${bigwig.version}/dist/BigWig.jar ${bigwig.jar}
	rm -rf bigwig-${bigwig.version}
	


##
## Apache Stuf
##
$(addprefix lib/, commons-validator/commons-validator/1.4.0/commons-validator-1.4.0.jar commons-beanutils/commons-beanutils/1.8.3/commons-beanutils-1.8.3.jar ):
	mkdir -p $(dir $@) && curl -Lk ${curl.proxy} -o $@ "http://central.maven.org/maven2/$(patsubst lib/%,%,$@)"

## API EVS
src/main/generated-sources/java/edu/washington/gs/evs/package-info.java : ${xjc_jars}
	mkdir -p ${generated.dir}/java
	${XJC} ${xjc.proxy} -d ${generated.dir}/java \
		-p edu.washington.gs.evs \
		"http://evs.gs.washington.edu/wsEVS/EVSDataQueryService?wsdl"

## copy resources for jxngs
jfxngs-resources:
	mkdir -p ${tmp.dir}/com/github/lindenb/jvarkit/tools/vcfviewgui
	$(foreach S,bam vcf, cp -v src/main/java/com/github/lindenb/jvarkit/tools/vcfviewgui/${S}.snippets.xml ${tmp.dir}/com/github/lindenb/jvarkit/tools/vcfviewgui/;)

#
# ga4gh schema, see 
# http://plindenbaum.blogspot.fr/2015/06/playing-with-ga4gh-schemas-and-avro-my.html
#
ga4gh.schemas.version ?= 0.5.1
ga4gh.schemas.avpr =  $(addsuffix .avpr,$(addprefix ${generated.dir}/resources/avro/ga4gh/schemas-${ga4gh.schemas.version}/src/main/resources/avro/,beacon common readmethods reads referencemethods references variantmethods variants wip/metadata wip/metadatamethods wip/variationReference ))

ga4gh.schemas.sources : lib/org/ga4gh/ga4gh-${ga4gh.schemas.version}.jar
lib/org/ga4gh/ga4gh-${ga4gh.schemas.version}.jar: ${avro.libs} ${ga4gh.schemas.avpr} 
	mkdir -p ${generated.dir}/java _tmp.avro${ga4gh.schemas.version} lib/org/ga4gh
	java -jar $< compile protocol \
		${generated.dir}/resources/avro/ga4gh/schemas-${ga4gh.schemas.version}/src/main/resources/avro \
		${generated.dir}/resources/avro/ga4gh/schemas-${ga4gh.schemas.version}/src/main/resources/avro/wip \
		${generated.dir}/java
	javac -d _tmp.avro${ga4gh.schemas.version} -cp ${avro.libs} -sourcepath ${generated.dir}/java ${generated.dir}/java/org/ga4gh/*.java
	jar cvf lib/org/ga4gh/ga4gh-${ga4gh.schemas.version}.jar -C _tmp.avro${ga4gh.schemas.version} .
	rm -rf _tmp.avro${ga4gh.schemas.version}
	 
	
$(filter-out $(word 1,${ga4gh.schemas.avpr}),${ga4gh.schemas.avpr}) : $(word 1,${ga4gh.schemas.avpr})
$(word 1,${ga4gh.schemas.avpr}) : ${avro.libs}
	rm -rf ${generated.dir}/resources/avro/ga4gh/schemas-${ga4gh.schemas.version}
	mkdir -p $(dir $@) 
	curl -Lk ${curl.proxy} -o schema.zip "https://github.com/ga4gh/schemas/archive/v${ga4gh.schemas.version}.zip"
	unzip schema.zip -d ${generated.dir}/resources/avro/ga4gh
	rm schema.zip	
	$(foreach P,${ga4gh.schemas.avpr}, java -jar ${avro.libs}  idl $(patsubst %.avpr,%.avdl,${P}) ${P} ; )
	touch $@

#
# end of ga4gh schema
#
#

#
# begin jni+htslib
#
#
HTSLIB_VERSION=1.6

htslib.jar : ${dist.dir}/libhtslibjni.so

${dist.dir}/libhtslibjni.so : ${dist.dir}/htslib.jar
	touch -c $@
${dist.dir}/htslib.jar :  $(addsuffix .java,$(addprefix ${src.dir}/com/github/lindenb/jvarkit/htslib/,HtsFile KString Bcf1)) \
		${this.dir}src/main/cpp/htslib/libhts.so \
		${this.dir}src/main/cpp/htslibjni.c
	rm -rf "${tmp.dir}" "${dist.dir}/libhtslibjni.so"
	mkdir -p $(dir $@)  "${tmp.dir}"
	$(JAVAC) -d ${tmp.dir}  -g -classpath "$(subst $(SPACE),:,$(filter %.jar,$^))" -sourcepath ${src.dir}:${generated.dir}/java $(filter %.java,$^)
	$(JAVAH) -o "${generated.dir}/cpp/htslibjni.h" -jni -classpath "${tmp.dir}" $(subst /,.,$(subst ${src.dir}/,,$(basename $(filter %.java,$^))))
	gcc -Wall -fPIC -shared -I /home/lindenb/package/jdk1.8.0_152/include -I /home/lindenb/package/jdk1.8.0_152/include/linux $(addprefix -I,$(sort $(dir $(shell find "${JAVA_HOME}/include" -type f -name "*.h")))) -I ${this.dir}src/main/cpp/htslib -I "${generated.dir}/cpp"  -o ${dist.dir}/libhtslibjni.so  $(filter %.c,$^)
	${JAR} cf $@ -C ${tmp.dir} .
	rm -rf "${tmp.dir}"
	

${this.dir}src/main/cpp/htslib/libhts.so :  ${this.dir}src/main/cpp/htslib/Makefile
	cd $(dir $@) && make
${this.dir}src/main/cpp/htslib/Makefile:
	rm -rf "$(dir $@)"
	wget -O "${HTSLIB_VERSION}.tar.gz" "https://github.com/samtools/htslib/archive/${HTSLIB_VERSION}.tar.gz"
	tar xfz "${HTSLIB_VERSION}.tar.gz"
	mv -v "htslib-${HTSLIB_VERSION}" "$(dir $@)"
	rm -vf "${HTSLIB_VERSION}.tar.gz"
	touch -c "$@"

#
# end of jni+htslib
#
#



include jfx.mk

## Knime helper
knimehelper: ${dist.dir}/knimehelper.jar
	@echo -n "KnimeHelper-";date +%Y-%m-%d
	@echo -n "KnimeHelper: A  java library to be used in the java nodes of [http://knime.org](http://knime.org). See [http://lindenb.github.io/jvarkit/KnimeIntegration.html](http://lindenb.github.io/jvarkit/KnimeIntegration.html). This library was compiled on "; LANG=en_UTF-8 date | tr -d "\n"; echo "."

${dist.dir}/knimehelper.jar: ${src.dir}/com/github/lindenb/jvarkit/knime/KnimeVariantHelper.java ${htsjdk.jars} ${jcommander.jar}
	mkdir -p ${tmp.dir}/META-INF/services
	${JAVAC} -cp "$(subst $(SPACE),:,$(realpath $(filter %.jar,$^)))" -d ${tmp.dir} -sourcepath ${src.dir}:${generated.dir}/java $<
	$(foreach J,$(filter %.jar,$^),unzip -o ${J} -d ${tmp.dir};) 
	${JAR} cvf $@ -C ${tmp.dir} .
	rm -rf "${tmp.dir}"


clean:
	rm -rf ${dist.dir}

