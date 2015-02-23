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
JAVAC?=javac
JAVA?=java
JAR?=jar
XJC?=xjc

htsjdk.version?=1.128
htsjdk.home?=${this.dir}htsjdk-${htsjdk.version}
htsjdk.jars=$(addprefix ${htsjdk.home}/dist/,$(addsuffix .jar,commons-jexl-2.1.1 commons-logging-1.1.1 htsjdk-${htsjdk.version} snappy-java-1.0.3-rc3))
src.dir=${this.dir}src/main/java
generated.dir=${this.dir}src/main/generated-sources
tmp.dir=${this.dir}_tmp-${htsjdk.version}
tmp.mft=${tmp.dir}/META-INF/MANIFEST.MF
dist.dir?=${this.dir}dist-${htsjdk.version}

mysql.version?=5.1.34
mysql.jar?=lib/mysql-connector-java-${mysql.version}-bin.jar
berkeleydb.version?=6.2.31
berkeleydb.jar?=lib/je-${berkeleydb.version}.jar
bigwig.jar?=lib/BigWig.jar


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
	echo "### COMPILING $(1) ######"
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


##
## Creates a CGI script for java app
##
define compile-cgi-cmd

${dist.dir}/$(1).cgi : $(1)
	echo "#!/bin/bash" > $$@
	echo "PREFIX=$$(dirname $$0)" >> $$@
	echo "${JAVA} -Xmx500m -Dprefs.file.xml=/var/www/cgi-bin/prefs.xml -jar $PREFIX/$(1).jar $$*" >> $$@

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
APPS= 		addlinearindextobed	allelefreqcalc	almostsortedvcf	backlocate	bam2fastq	bam2raster	bam2svg \
	bam2wig	bam4deseq01	bamcmpcoverage	bamgenscan	bamindexreadnames	bamliftover	bamqueryreadnames \
	bamrenamechr	bamsnvwig	bamstats04	bamtreepack	bamviewgui	batchigvpictures	bedliftover \
	bedrenamechr	biostar103303	biostar106668	biostar130456	biostar59647	biostar76892	biostar77288 \
	biostar77828	biostar78285	biostar78400	biostar81455	biostar84452	biostar84786	biostar86363 \
	biostar86480	biostar90204	biostar94573	95652	blast2sam	blastfastq	blastmapannots \
	blastn2snp	buildwpontology	bwamemdigest	bwamemnop	cmpbams	cmpbamsandbuild	coveragenormalizer \
	deseqcount	downsamplevcf	evs2bed	evs2vcf	evs2xml	extendbed	fastq2fasta \
	fastqentropy	fastqgrep	fastqjs	fastqphred64to33	fastqrecordtreepack	fastqrevcomp	fastqshuffle \
	fastqsplitinterleaved	findallcoverageatposition	findavariation	findcorruptedfiles	findmyvirus	findnewsplicesites	fixvarscanmissingheader \
	fixvcf	fixvcfformat	fixvcfmissinggenotypes	gcanddepth	genomicjaspar	genscan	groupbygene \
	howmanybamdict	idea20130924	illuminadir	ilmnfastqstats	impactofduplicates	jeter	kg2bed \
	liftover2svg	mapuniprot	mergesplittedblast	metrics2xml	ncbitaxonomy2xml	ngsfilessummary	noemptyvcf \
	nozerovariationvcf	pademptyfastq	paintcontext	pubmeddump	pubmedfilterjs	referencetovcf	sam2json \
	sam2psl	sam2tsv	sam4weblogo	samclipindelfraction	samextractclip	samfindclippedregions	samfixcigar \
	samgrep	samjs	samshortinvert	samstats01	scanshortinvert	sigframe	sortvcfoninfo \
	sortvcfonref2	splitbam	splitbam2	splitbytile	splitread	tview	tview.cgi \
	vcf2hilbert	vcf2ps	vcf2rdf	vcf2sql	vcf2xml	vcfannobam	vcfbed \
	vcfbedjs	vcfbiomart	vcfcadd	vcfcmppred	vcfcomm	vcfcompare	vcfcomparegt \
	vcfconcat	vcfcutsamples	vcffilterdoid	vcffilterjs	vcffilterso	vcffixindels	vcfgo \
	vcfhead	vcfin	vcfjaspar	vcfliftover	vcfmapuniprot	vcfmerge	vcfmulti2one \
	vcfpolyx	vcfpredictions	vcfrebase	vcfregistry.cgi	vcfregulomedb	vcfrenamechr	vcfrenamesamples \
	vcfresetvcf	vcfsetdict	vcfshuffle	vcfsimulator	vcfstats	vcfstopcodon	vcfstripannot \
	vcftabixml	vcftail	vcftreepack	vcftrio	vcfvcf	vcfviewgui	worldmapgenome \


.PHONY: all $(APPS) clean library

all: $(APPS)



ifneq ($(realpath ${bigwig.jar}),)
$(eval $(call compile-htsjdk-cmd,vcfbigwig,		com.github.lindenb.jvarkit.tools.vcfbigwig.VCFBigWig,${bigwig.jar}))
$(eval $(call compile-htsjdk-cmd,vcfensemblreg,	com.github.lindenb.jvarkit.tools.ensemblreg.VcfEnsemblReg,${bigwig.jar}))
$(eval $(call compile_biostar_cmd,105754,${bigwig.jar}))
endif

ifneq ($(realpath ${common.math.jar}),)
$(eval $(call compile-htsjdk-cmd,cnv01,com.github.lindenb.jvarkit.tools.redon.CopyNumber01,${common.math.jar}))
endif

ifneq ($(realpath ${berkeleydb.jar}),)
$(eval $(call compile-htsjdk-cmd,vcfphylotree,com.github.lindenb.jvarkit.tools.phylo.VcfPhyloTree,${berkeleydb.jar}))
$(eval $(call compile-htsjdk-cmd,ngsfilesscanner,com.github.lindenb.jvarkit.tools.ngsfiles.NgsFilesScanner,${berkeleydb.jar}))
$(eval $(call compile-htsjdk-cmd,92368,${berkeleydb.jar}))
endif

ifneq ($(realpath ${mysql.jar}),)
$(eval $(call compile-htsjdk-cmd,vcfucsc,com.github.lindenb.jvarkit.tools.vcfucsc.VcfUcsc,${mysql.jar}))
endif

$(eval $(call compile-htsjdk-cmd,addlinearindextobed,com.github.lindenb.jvarkit.tools.misc.AddLinearIndexToBed))
$(eval $(call compile-htsjdk-cmd,allelefreqcalc,com.github.lindenb.jvarkit.tools.misc.AlleleFrequencyCalculator))
$(eval $(call compile-htsjdk-cmd,almostsortedvcf,com.github.lindenb.jvarkit.tools.sortvcfonref.AlmostSortedVcf))
$(eval $(call compile-htsjdk-cmd,backlocate,com.github.lindenb.jvarkit.tools.backlocate.BackLocate))
$(eval $(call compile-htsjdk-cmd,bam2fastq,com.github.lindenb.jvarkit.tools.fastq.BamToFastq))
$(eval $(call compile-htsjdk-cmd,bam2raster,com.github.lindenb.jvarkit.tools.bam2graphics.Bam2Raster))
$(eval $(call compile-htsjdk-cmd,bam2svg,com.github.lindenb.jvarkit.tools.bam2svg.BamToSVG))
$(eval $(call compile-htsjdk-cmd,bam2wig,com.github.lindenb.jvarkit.tools.bam2wig.Bam2Wig))
$(eval $(call compile-htsjdk-cmd,bam4deseq01,com.github.lindenb.jvarkit.tools.bam4deseq.Bam4DeseqIntervals))
$(eval $(call compile-htsjdk-cmd,bamcmpcoverage,com.github.lindenb.jvarkit.tools.misc.BamCmpCoverage))
$(eval $(call compile-htsjdk-cmd,bamgenscan,com.github.lindenb.jvarkit.tools.genscan.BamGenScan))
$(eval $(call compile-htsjdk-cmd,bamindexreadnames,com.github.lindenb.jvarkit.tools.bamindexnames.BamIndexReadNames))
$(eval $(call compile-htsjdk-cmd,bamliftover,com.github.lindenb.jvarkit.tools.liftover.BamLiftOver))
$(eval $(call compile-htsjdk-cmd,bamqueryreadnames,com.github.lindenb.jvarkit.tools.bamindexnames.BamQueryReadNames))
$(eval $(call compile-htsjdk-cmd,bamrenamechr,com.github.lindenb.jvarkit.tools.misc.ConvertBamChromosomes))
$(eval $(call compile-htsjdk-cmd,bamsnvwig,com.github.lindenb.jvarkit.tools.mem.BWAMemScan))
$(eval $(call compile-htsjdk-cmd,bamstats04,com.github.lindenb.jvarkit.tools.bamstats04.BamStats04))
$(eval $(call compile-htsjdk-cmd,bamtreepack,com.github.lindenb.jvarkit.tools.treepack.BamTreePack))
$(eval $(call compile-htsjdk-cmd,bamviewgui,com.github.lindenb.jvarkit.tools.bamviewgui.BamViewGui))
$(eval $(call compile-htsjdk-cmd,batchigvpictures,com.github.lindenb.jvarkit.tools.batchpicts.BatchIGVPictures,copy.opendoc.odp.resources))
$(eval $(call compile-htsjdk-cmd,bedliftover,com.github.lindenb.jvarkit.tools.liftover.BedLiftOver))
$(eval $(call compile-htsjdk-cmd,bedrenamechr,com.github.lindenb.jvarkit.tools.misc.ConvertBedChromosomes))
$(eval $(call compile_biostar_cmd,103303))
$(eval $(call compile_biostar_cmd,106668))
$(eval $(call compile_biostar_cmd,130456))
$(eval $(call compile_biostar_cmd,59647))
$(eval $(call compile_biostar_cmd,76892))
$(eval $(call compile_biostar_cmd,77288))
$(eval $(call compile_biostar_cmd,77828))
$(eval $(call compile_biostar_cmd,78285))
$(eval $(call compile_biostar_cmd,78400))
$(eval $(call compile_biostar_cmd,81455))
$(eval $(call compile_biostar_cmd,84452))
$(eval $(call compile_biostar_cmd,84786))
$(eval $(call compile_biostar_cmd,86363))
$(eval $(call compile_biostar_cmd,86480))
$(eval $(call compile_biostar_cmd,90204))
$(eval $(call compile_biostar_cmd,94573))
$(eval $(call compile_biostar_cmd,95652,api.ncbi.gb))
$(eval $(call compile-htsjdk-cmd,blast2sam,com.github.lindenb.jvarkit.tools.blast2sam.BlastToSam,api.ncbi.blast))
$(eval $(call compile-htsjdk-cmd,blastfastq,com.github.lindenb.jvarkit.tools.bwamempcr.BlastFastQ))
$(eval $(call compile-htsjdk-cmd,blastmapannots, com.github.lindenb.jvarkit.tools.blastmapannots.BlastMapAnnotations, api.ncbi.blast,api.ncbi.gb api.uniprot))
$(eval $(call compile-htsjdk-cmd,blastn2snp,com.github.lindenb.jvarkit.tools.blast.BlastNToSnp,api.ncbi.blast))
$(eval $(call compile-htsjdk-cmd,buildwpontology,com.github.lindenb.jvarkit.tools.misc.BuildWikipediaOntology))
$(eval $(call compile-htsjdk-cmd,bwamemdigest,com.github.lindenb.jvarkit.tools.mem.BWAMemDigest))
$(eval $(call compile-htsjdk-cmd,bwamemnop,com.github.lindenb.jvarkit.tools.mem.BWAMemNOp))
$(eval $(call compile-htsjdk-cmd,cmpbams,com.github.lindenb.jvarkit.tools.cmpbams.CompareBams2))
$(eval $(call compile-htsjdk-cmd,cmpbamsandbuild,com.github.lindenb.jvarkit.tools.cmpbams.CompareBamAndBuild))
$(eval $(call compile-htsjdk-cmd,coveragenormalizer,com.github.lindenb.jvarkit.tools.misc.CoverageNormalizer))
$(eval $(call compile-htsjdk-cmd,deseqcount,com.github.lindenb.jvarkit.tools.bam4deseq.HtSeqCount))
$(eval $(call compile-htsjdk-cmd,downsamplevcf,com.github.lindenb.jvarkit.tools.misc.DownSampleVcf))
$(eval $(call compile-htsjdk-cmd,evs2bed,com.github.lindenb.jvarkit.tools.evs2bed.DumpExomeVariantServerData))
$(eval $(call compile-htsjdk-cmd,evs2vcf,com.github.lindenb.jvarkit.tools.evs2bed.EvsToVcf,api.evs))
$(eval $(call compile-htsjdk-cmd,evs2xml,com.github.lindenb.jvarkit.tools.evs2bed.EvsDumpXml,api.evs))
$(eval $(call compile-htsjdk-cmd,extendbed,com.github.lindenb.jvarkit.tools.misc.ExtendBed))
$(eval $(call compile-htsjdk-cmd,fastq2fasta,com.github.lindenb.jvarkit.tools.misc.FastqToFasta))
$(eval $(call compile-htsjdk-cmd,fastqentropy,com.github.lindenb.jvarkit.tools.fastq.FastqEntropy))
$(eval $(call compile-htsjdk-cmd,fastqgrep,com.github.lindenb.jvarkit.tools.misc.FastqGrep))
$(eval $(call compile-htsjdk-cmd,fastqjs,com.github.lindenb.jvarkit.tools.fastq.FastqJavascript))
$(eval $(call compile-htsjdk-cmd,fastqphred64to33,com.github.lindenb.jvarkit.tools.fastq.ConvertPhred64toFastq33))
$(eval $(call compile-htsjdk-cmd,fastqrecordtreepack,com.github.lindenb.jvarkit.tools.treepack.FastqRecordTreePack))
$(eval $(call compile-htsjdk-cmd,fastqrevcomp,com.github.lindenb.jvarkit.tools.misc.FastqRevComp))
$(eval $(call compile-htsjdk-cmd,fastqshuffle,com.github.lindenb.jvarkit.tools.fastq.FastqShuffle))
$(eval $(call compile-htsjdk-cmd,fastqsplitinterleaved,com.github.lindenb.jvarkit.tools.fastq.FastqSplitInterleaved))
$(eval $(call compile-htsjdk-cmd,findallcoverageatposition,com.github.lindenb.jvarkit.tools.misc.FindAllCoverageAtPosition))
$(eval $(call compile-htsjdk-cmd,findavariation,com.github.lindenb.jvarkit.tools.misc.FindAVariation))
$(eval $(call compile-htsjdk-cmd,findcorruptedfiles,com.github.lindenb.jvarkit.tools.misc.FindCorruptedFiles))
$(eval $(call compile-htsjdk-cmd,findmyvirus,com.github.lindenb.jvarkit.tools.mem.FindMyVirus))
$(eval $(call compile-htsjdk-cmd,findnewsplicesites,com.github.lindenb.jvarkit.tools.rnaseq.FindNewSpliceSites))
$(eval $(call compile-htsjdk-cmd,fixvarscanmissingheader,com.github.lindenb.jvarkit.tools.misc.FixVarScanMissingVCFHeader))
$(eval $(call compile-htsjdk-cmd,fixvcf,com.github.lindenb.jvarkit.tools.misc.FixVCF))
$(eval $(call compile-htsjdk-cmd,fixvcfformat,com.github.lindenb.jvarkit.tools.misc.FixVcfFormat))
$(eval $(call compile-htsjdk-cmd,fixvcfmissinggenotypes,com.github.lindenb.jvarkit.tools.misc.FixVcfMissingGenotypes))
$(eval $(call compile-htsjdk-cmd,gcanddepth,com.github.lindenb.jvarkit.tools.misc.GcPercentAndDepth))
$(eval $(call compile-htsjdk-cmd,genomicjaspar,com.github.lindenb.jvarkit.tools.jaspar.GenomicJaspar))
$(eval $(call compile-htsjdk-cmd,genscan,com.github.lindenb.jvarkit.tools.genscan.GenScan))
$(eval $(call compile-htsjdk-cmd,groupbygene,com.github.lindenb.jvarkit.tools.groupbygene.GroupByGene))
$(eval $(call compile-htsjdk-cmd,howmanybamdict,com.github.lindenb.jvarkit.tools.misc.HowManyBamDict))
$(eval $(call compile-htsjdk-cmd,idea20130924,com.github.lindenb.jvarkit.tools.bwamempcr.Idea20130924))
$(eval $(call compile-htsjdk-cmd,illuminadir,com.github.lindenb.jvarkit.tools.misc.IlluminaDirectory))
$(eval $(call compile-htsjdk-cmd,ilmnfastqstats,com.github.lindenb.jvarkit.tools.misc.IlluminaStatsFastq))
$(eval $(call compile-htsjdk-cmd,impactofduplicates,com.github.lindenb.jvarkit.tools.impactdup.ImpactOfDuplicates))
$(eval $(call compile-htsjdk-cmd,kg2bed,com.github.lindenb.jvarkit.tools.misc.KnownGenesToBed))
$(eval $(call compile-htsjdk-cmd,liftover2svg,com.github.lindenb.jvarkit.tools.liftover.LiftOverToSVG))
$(eval $(call compile-htsjdk-cmd,mapuniprot,com.github.lindenb.jvarkit.tools.misc.MapUniProtFeatures,api.uniprot))
$(eval $(call compile-htsjdk-cmd,mergesplittedblast,com.github.lindenb.jvarkit.tools.blast.MergeSplittedBlast,api.ncbi.blast))
$(eval $(call compile-htsjdk-cmd,metrics2xml,com.github.lindenb.jvarkit.tools.metrics2xml.PicardMetricsToXML))
$(eval $(call compile-htsjdk-cmd,ncbitaxonomy2xml,com.github.lindenb.jvarkit.tools.misc.NcbiTaxonomyToXml))
$(eval $(call compile-htsjdk-cmd,ngsfilessummary,com.github.lindenb.jvarkit.tools.ngsfiles.NgsFilesSummary))
$(eval $(call compile-htsjdk-cmd,noemptyvcf,com.github.lindenb.jvarkit.tools.misc.NoEmptyVCF))
$(eval $(call compile-htsjdk-cmd,nozerovariationvcf,com.github.lindenb.jvarkit.tools.misc.NoZeroVariationVCF))
$(eval $(call compile-htsjdk-cmd,pademptyfastq,com.github.lindenb.jvarkit.tools.misc.PadEmptyFastq))
$(eval $(call compile-htsjdk-cmd,paintcontext,com.github.lindenb.jvarkit.tools.bam2graphics.PaintContext))
$(eval $(call compile-htsjdk-cmd,pubmeddump,com.github.lindenb.jvarkit.tools.pubmed.PubmedDump))
$(eval $(call compile-htsjdk-cmd,pubmedfilterjs,com.github.lindenb.jvarkit.tools.pubmed.PubmedFilterJS,api.ncbi.pubmed))
$(eval $(call compile-htsjdk-cmd,referencetovcf,com.github.lindenb.jvarkit.tools.misc.ReferenceToVCF))
$(eval $(call compile-htsjdk-cmd,sam2json,com.github.lindenb.jvarkit.tools.misc.SamToJson))
$(eval $(call compile-htsjdk-cmd,sam2psl,com.github.lindenb.jvarkit.tools.misc.SamToPsl))
$(eval $(call compile-htsjdk-cmd,sam2tsv,com.github.lindenb.jvarkit.tools.sam2tsv.Sam2Tsv))
$(eval $(call compile-htsjdk-cmd,sam4weblogo,com.github.lindenb.jvarkit.tools.sam4weblogo.SAM4WebLogo))
$(eval $(call compile-htsjdk-cmd,samclipindelfraction,com.github.lindenb.jvarkit.tools.misc.SamClipIndelFraction))
$(eval $(call compile-htsjdk-cmd,samextractclip,com.github.lindenb.jvarkit.tools.structvar.SamExtractClip))
$(eval $(call compile-htsjdk-cmd,samfindclippedregions,com.github.lindenb.jvarkit.tools.structvar.SamFindClippedRegions))
$(eval $(call compile-htsjdk-cmd,samfixcigar,com.github.lindenb.jvarkit.tools.samfixcigar.SamFixCigar))
$(eval $(call compile-htsjdk-cmd,samgrep,com.github.lindenb.jvarkit.tools.samgrep.SamGrep))
$(eval $(call compile-htsjdk-cmd,samjs,com.github.lindenb.jvarkit.tools.samjs.SamJavascript))
$(eval $(call compile-htsjdk-cmd,samshortinvert,com.github.lindenb.jvarkit.tools.structvar.SamShortInvertion))
$(eval $(call compile-htsjdk-cmd,samstats01,com.github.lindenb.jvarkit.tools.bamstats01.BamStats01))
$(eval $(call compile-htsjdk-cmd,scanshortinvert,com.github.lindenb.jvarkit.tools.mem.ScanShortInvert))
$(eval $(call compile-htsjdk-cmd,sigframe,com.github.lindenb.jvarkit.tools.sigframe.SigFrame))
$(eval $(call compile-htsjdk-cmd,sortvcfoninfo,com.github.lindenb.jvarkit.tools.sortvcfonref.SortVcfOnInfo))
$(eval $(call compile-htsjdk-cmd,sortvcfonref2,com.github.lindenb.jvarkit.tools.sortvcfonref.SortVcfOnRef2))
$(eval $(call compile-htsjdk-cmd,splitbam,com.github.lindenb.jvarkit.tools.splitbam.SplitBam))
$(eval $(call compile-htsjdk-cmd,splitbam2,com.github.lindenb.jvarkit.tools.splitbam.SplitBam2))
$(eval $(call compile-htsjdk-cmd,splitbytile,com.github.lindenb.jvarkit.tools.splitbytitle.SplitByTile))
$(eval $(call compile-htsjdk-cmd,splitread,com.github.lindenb.jvarkit.tools.splitread.SplitRead))
#$(eval $(call compile-htsjdk-cmd,tview,com.github.lindenb.jvarkit.tools.tview.TViewCmd))
#$(eval $(call compile-cgi-cmd,tview.cgi))
$(eval $(call compile-htsjdk-cmd,vcf2hilbert,com.github.lindenb.jvarkit.tools.misc.VcfToHilbert))
$(eval $(call compile-htsjdk-cmd,vcf2ps,com.github.lindenb.jvarkit.tools.misc.VcfToPostscript))
$(eval $(call compile-htsjdk-cmd,vcf2rdf,com.github.lindenb.jvarkit.tools.vcf2rdf.VcfToRdf))
$(eval $(call compile-htsjdk-cmd,vcf2sql,com.github.lindenb.jvarkit.tools.vcf2sql.VcfToSql))
$(eval $(call compile-htsjdk-cmd,vcf2xml,com.github.lindenb.jvarkit.tools.vcf2xml.Vcf2Xml))
$(eval $(call compile-htsjdk-cmd,vcfannobam,com.github.lindenb.jvarkit.tools.vcfannobam.VCFAnnoBam))
$(eval $(call compile-htsjdk-cmd,vcfbed,com.github.lindenb.jvarkit.tools.vcfbed.VCFBed))
$(eval $(call compile-htsjdk-cmd,vcfbedjs,com.github.lindenb.jvarkit.tools.vcfbed.VCFBed))
$(eval $(call compile-htsjdk-cmd,vcfbiomart,com.github.lindenb.jvarkit.tools.vcfbiomart.VcfBiomart))
$(eval $(call compile-htsjdk-cmd,vcfcadd,com.github.lindenb.jvarkit.tools.misc.VcfCadd))
$(eval $(call compile-htsjdk-cmd,vcfcmppred,com.github.lindenb.jvarkit.tools.vcfcmp.VCFComparePredictions))
$(eval $(call compile-htsjdk-cmd,vcfcomm,com.github.lindenb.jvarkit.tools.vcfcmp.VCFComm))
$(eval $(call compile-htsjdk-cmd,vcfcompare,com.github.lindenb.jvarkit.tools.vcfcmp.VCFCompare))
$(eval $(call compile-htsjdk-cmd,vcfcomparegt,com.github.lindenb.jvarkit.tools.vcfcmp.VCFCompareGT))
$(eval $(call compile-htsjdk-cmd,vcfconcat,com.github.lindenb.jvarkit.tools.vcfconcat.VcfConcat))
$(eval $(call compile-htsjdk-cmd,vcfcutsamples,com.github.lindenb.jvarkit.tools.misc.VcfCutSamples))
$(eval $(call compile-htsjdk-cmd,vcfdas,com.github.lindenb.jvarkit.tools.vcfdas.VcfDistributedAnnotationSystem, ${jetty.jars}))
$(eval $(call compile-htsjdk-cmd,vcffilterdoid,com.github.lindenb.jvarkit.tools.vcfdo.VcfFilterDoid))
$(eval $(call compile-htsjdk-cmd,vcffilterjs,com.github.lindenb.jvarkit.tools.vcffilterjs.VCFFilterJS))
$(eval $(call compile-htsjdk-cmd,vcffilterso,com.github.lindenb.jvarkit.tools.misc.VcfFilterSequenceOntology))
$(eval $(call compile-htsjdk-cmd,vcffixindels,com.github.lindenb.jvarkit.tools.vcffixindels.VCFFixIndels))
$(eval $(call compile-htsjdk-cmd,vcfgo,com.github.lindenb.jvarkit.tools.vcfgo.VcfGeneOntology))
$(eval $(call compile-htsjdk-cmd,vcfhead,com.github.lindenb.jvarkit.tools.misc.VcfHead))
$(eval $(call compile-htsjdk-cmd,vcfin,com.github.lindenb.jvarkit.tools.vcfcmp.VcfIn))
$(eval $(call compile-htsjdk-cmd,vcfjaspar,com.github.lindenb.jvarkit.tools.jaspar.VcfJaspar))
$(eval $(call compile-htsjdk-cmd,vcfliftover,com.github.lindenb.jvarkit.tools.liftover.VcfLiftOver))
$(eval $(call compile-htsjdk-cmd,vcfmapuniprot,com.github.lindenb.jvarkit.tools.misc.VcfMapUniprot,api.uniprot))
$(eval $(call compile-htsjdk-cmd,vcfmerge,com.github.lindenb.jvarkit.tools.vcfmerge.VCFMerge2))
$(eval $(call compile-htsjdk-cmd,vcfmulti2one,com.github.lindenb.jvarkit.tools.misc.VcfMultiToOne))
$(eval $(call compile-htsjdk-cmd,vcfpolyx,com.github.lindenb.jvarkit.tools.misc.VCFPolyX))
$(eval $(call compile-htsjdk-cmd,vcfpredictions,com.github.lindenb.jvarkit.tools.vcfannot.VCFAnnotator))
$(eval $(call compile-htsjdk-cmd,vcfrebase,com.github.lindenb.jvarkit.tools.vcfrebase.VcfRebase))
$(eval $(call compile-cgi-cmd,vcfregistry.cgi))
$(eval $(call compile-htsjdk-cmd,vcfregulomedb,com.github.lindenb.jvarkit.tools.misc.VcfRegulomeDB))
$(eval $(call compile-htsjdk-cmd,vcfrenamechr,com.github.lindenb.jvarkit.tools.misc.ConvertVcfChromosomes))
$(eval $(call compile-htsjdk-cmd,vcfrenamesamples,com.github.lindenb.jvarkit.tools.misc.VcfRenameSamples))
$(eval $(call compile-htsjdk-cmd,vcfresetvcf,com.github.lindenb.jvarkit.tools.misc.VcfRemoveGenotypeIfInVcf))
$(eval $(call compile-htsjdk-cmd,vcfsetdict,com.github.lindenb.jvarkit.tools.misc.VcfSetSequenceDictionary))
$(eval $(call compile-htsjdk-cmd,vcfshuffle,com.github.lindenb.jvarkit.tools.misc.VCFShuffle))
$(eval $(call compile-htsjdk-cmd,vcfsimulator,com.github.lindenb.jvarkit.tools.misc.VcfSimulator))
$(eval $(call compile-htsjdk-cmd,vcfstats,com.github.lindenb.jvarkit.tools.vcfstats.VcfStats))
$(eval $(call compile-htsjdk-cmd,vcfstopcodon,com.github.lindenb.jvarkit.tools.vcfannot.VCFStopCodon))
$(eval $(call compile-htsjdk-cmd,vcfstripannot,com.github.lindenb.jvarkit.tools.vcfstripannot.VCFStripAnnotations))
$(eval $(call compile-htsjdk-cmd,vcftabixml,com.github.lindenb.jvarkit.tools.vcftabixml.VCFTabixml))
$(eval $(call compile-htsjdk-cmd,vcftail,com.github.lindenb.jvarkit.tools.misc.VcfTail))
$(eval $(call compile-htsjdk-cmd,vcftreepack,com.github.lindenb.jvarkit.tools.treepack.VcfTreePack))
$(eval $(call compile-htsjdk-cmd,vcftrio,com.github.lindenb.jvarkit.tools.vcftrios.VCFTrios))
$(eval $(call compile-htsjdk-cmd,vcfvcf,com.github.lindenb.jvarkit.tools.vcfvcf.VcfVcf))
$(eval $(call compile-htsjdk-cmd,vcfviewgui,com.github.lindenb.jvarkit.tools.vcfviewgui.VcfViewGui))
$(eval $(call compile-htsjdk-cmd,worldmapgenome,com.github.lindenb.jvarkit.tools.circular.WorldMapGenome))




all-jnlp : $(addprefix ${dist.dir}/,$(addsuffix .jar,vcfviewgui buildwpontology batchigvpictures)) ${htsjdk.jars} \
	 ./src/main/resources/jnlp/generic.jnlp .secret.keystore 
	mkdir -p ${tmp.dir}
	cp $(filter %.jar,$^) ${tmp.dir}
	$(foreach J, $(filter %.jar,$^) , jarsigner ${J} secret ; ) 
	sed -e 's/__MAIN_JAR__/vcfviewgui/g' \
		-e 's/__TITLE__/VcfViewGui/g' \
		-e 's/__MAIN_CLASS__/com.github.lindenb.jvarkit.tools.vcfviewgui.VcfViewGui/g' \
		 ./src/main/resources/jnlp/generic.jnlp > ${tmp.dir}/vcfviewgui.jnlp
	sed -e 's/__MAIN_JAR__/buildwpontology/g' \
		-e 's/__TITLE__/BuildWikipediaOntology/g' \
		-e 's/__MAIN_CLASS__/com.github.lindenb.jvarkit.tools.misc.BuildWikipediaOntology/g' \
		 ./src/main/resources/jnlp/generic.jnlp > ${tmp.dir}/buildwpontology.jnlp
	sed -e 's/__MAIN_JAR__/batchigvpictures/g' \
		-e 's/__TITLE__/BatchIGVPictures/g' \
		-e 's/__MAIN_CLASS__/com.github.lindenb.jvarkit.tools.batchpicts.BatchIGVPictures/g' \
		 ./src/main/resources/jnlp/generic.jnlp > ${tmp.dir}/batchigvpictures.jnlp

.secret.keystore : 
	 -keytool -genkeypair -keystore $@ -alias secret \
	 		-keypass $(if ${keytool.keypass},${keytool.keypass},KEYTOOLPASS) \
	 		-storepass $(if ${keytool.storepass},${keytool.storepass},KEYTOOLSTOREPASS) \
	 		-dname CN=Pierre Lindenbaum, OU=INSERM, O=INSERM, L=Nantes, ST=Nantes, C=Fr

api.uniprot :
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java -p org.uniprot ${xjc.proxy} "http://www.uniprot.org/support/docs/uniprot.xsd"


api.ncbi.pubmed : 
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java  -p gov.nih.nlm.ncbi.pubmed -dtd ${xjc.proxy} http://www.ncbi.nlm.nih.gov/corehtml/query/DTD/pubmed_140101.dtd


api.evs:
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java -p edu.washington.gs.evs ${xjc.proxy} -wsdl "http://evs.gs.washington.edu/wsEVS/EVSDataQueryService?wsdl"

api.ncbi.blast:
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java  -p gov.nih.nlm.ncbi.blast -dtd ${xjc.proxy} "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd"

api.ncbi.esearch:
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java  -p gov.nih.nlm.ncbi.esearch -dtd ${xjc.proxy} "http://eutils.ncbi.nlm.nih.gov/eutils/dtd/20060628/esearch.dtd"

api.ncbi.gb:
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java  -p gov.nih.nlm.ncbi.esearch -dtd ${xjc.proxy} "http://www.ncbi.nlm.nih.gov/dtd/NCBI_GBSeq.dtd"

api.ncbi.taxonomy:
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java  -p gov.nih.nlm.ncbi.taxonomy -dtd ${xjc.proxy} "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/taxon.dtd"

api.ncbi.tseq:
	mkdir -p ${generated.dir}/java
	${XJC} -d ${generated.dir}/java  -p gov.nih.nlm.ncbi.tseq -dtd ${xjc.proxy} "http://www.ncbi.nlm.nih.gov/dtd/NCBI_TSeq.dtd"
	
	
api.printf:
	mkdir -p ${generated.dir}/java/com/github/lindenb/jvarkit/util/printf
	javacc -OUTPUT_DIRECTORY=${generated.dir}/java/com/github/lindenb/jvarkit/util/printf \
		src/main/resources/javacc/com/github/lindenb/jvarkit/util/printf/Printf.jj

copy.opendoc.odp.resources : 
	mkdir -p ${tmp.dir}/META-INF/opendocument/odp
	cp src/main/resources/opendocument/thumbnail.png ${tmp.dir}/META-INF/opendocument/odp/
	cp $(addprefix src/main/resources/opendocument/odp/, meta.xml content.xml settings.xml styles.xml mimetype ) ${tmp.dir}/META-INF/opendocument/odp/


## jvarkit-library (used in knime)
library: ${dist.dir}/jvarkit-${htsjdk.version}.jar
${dist.dir}/jvarkit-${htsjdk.version}.jar : ${htsjdk.jars} \
		${generated.dir}/java/com/github/lindenb/jvarkit/util/htsjdk/HtsjdkVersion.java \
		${src.dir}/com/github/lindenb/jvarkit/util/Library.java
	mkdir -p ${tmp.dir}/META-INF $(dir $@)
	cp src/main/resources/messages/messages.properties ${tmp.dir}
	${JAVAC} -d ${tmp.dir} -g -classpath "$(subst $(SPACE),:,$(filter %.jar,$^))" -sourcepath ${src.dir}:${generated.dir}/java $(filter %.java,$^)
	${JAR} cf $@ -C ${tmp.dir} .
	rm -rf ${tmp.dir}

##
## Download mysql connector for java
## 
${mysql.jar} :
	echo "Downloading mysql connector for java version ${mysql.version} from oracle"
	mkdir -p $(dir $@)
	curl -Lk ${curl.proxy} -o jeter.tar.gz "http://dev.mysql.com/get/Downloads/Connector-J/mysql-connector-java-${mysql.version}.tar.gz" 
	tar xvfz  jeter.tar.gz  --strip-components=1 "mysql-connector-java-${mysql.version}/mysql-connector-java-${mysql.version}-bin.jar"
	mv mysql-connector-java-${mysql.version}-bin.jar $@
	rm jeter.tar.gz

##
## ORACLE BerkeleyDB for java
## 
${berkeleydb.jar}:
	echo "Downloading berkeleydb for java version ${berkeleydb.version} from oracle"
	mkdir -p $($@)
	curl -Lk ${curl.proxy} -o $@ "http://download.oracle.com/maven/com/sleepycat/je/${berkeleydb.version}/je-${berkeleydb.version}.jar"

${bigwig.jar} :
	echo "Downloading bigwig library for java. Requires subversion (SVN) and apache ANT"
	mkdir -p $(dir $@)
	rm -rf bigwig-read-only 
	svn checkout "http://bigwig.googlecode.com/svn/trunk/" bigwig-read-only
	(cd bigwig-read-only; ant)
	mv bigwig-read-only/dist/BigWig.jar $@
	rm -rf bigwig-read-only 
	

##
## make sure jars from htslib exist
##
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

test-vcfgo: vcfgo
	echo "" | $(JAVA) -jar ${dist.dir}/vcfgo.jar \
		-A "http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.goa_human.gz?rev=HEAD" \
		-G "http://archive.geneontology.org/latest-termdb/go_daily-termdb.rdf-xml.gz" \
		 -C GO:0007283 -F GOFILTER -v


clean:
	rm -rf ${dist.dir}


