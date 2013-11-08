JVARKIT
=======

Java utilities for Next Generation Sequencing

Pierre Lindenbaum PhD

http://plindenbaum.blogspot.com

@yokofakun	

## Download and install

see [Download and Install](https://github.com/lindenb/jvarkit/wiki/Compilation)

<h2>Tools</h2>
I'm slowly moving the pages below to the github wiki  ( https://github.com/lindenb/jvarkit/wiki/_pages )
<table>
<tr><th>Tool</th><th>Description</th></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Illuminadir">Illuminadir</a></th><td>Create a structured (**JSON** or **XML**) representation of a directory containing some Illumina 
FASTQs.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SamJS">SamJS</a></th><td>Filtering a SAM/BAM with javascript (rhino).</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFFilterJS">VCFFilterJS</a></th><td>Filtering a VCF with javascript (rhino)</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SplitBam">SplitBam</a></th><td>Split a BAM by chromosome group. Creates EMPTY bams if no reads was found for a given group. </td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCF2XML">VCF2XML</a></th><td>Transforms a VCF to XML. </td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFPredictions">VCFPredictions<a></th><td>Basic variant prediction using UCSC knownGenes.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SortVCFOnRef">SortVCFOnRef<a></th><td>Sort a VCF using the order of the chromosomes in a REFerence index.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SAM4WebLogo">SAM4WebLogo<a></th><td>Creates an Input file for BAM + WebLogo.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SAM2Tsv">SAM2Tsv<a></th><td>Tabular view of each base of the reads vs the reference.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/CmpBams">CmpBams<a></th><td>Compare two or more BAMs.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar84786">Biostar84786<a></th><td>Table transposition</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/IlluminaFastqStats">IlluminaFastqStats<a></th><td>Statistics on Illumina Fastqs</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Bam2Raster">Bam2Raster<a></th><td>Save a BAM alignment as a PNG image.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCF2SQL">VCF2SQL<a></th><td>Generate the SQL code to insert a VCF into a database</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Bam4DeseqIntervals">Bam4DeseqIntervals<a></th><td>creates a table for DESEQ with the number of reads within a sliding window for multiple BAMS</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFStripAnnotations">VCFStripAnnotations<a></th><td>Removes one or more field from the INFO column from a VCF.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFFixIndels">VCFFixIndels<a></th><td>Fix samtools INDELS for @SolenaLS</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BlastMapAnnots">BlastMapAnnots<a></th><td>Maps uniprot/genbank annotations on a blast result. See http://www.biostars.org/p/76056</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfViewGui">VcfViewGui<a></th><td>Simple java-Swing-based VCF viewer.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SamGrep">SamGrep<a></th><td>Search reads in a BAM</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFGeneOntology">VCFGeneOntology<a></th><td>Finds the GO terms for VCF annotated with SNPEFF or VEP</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BamStats01">BamStats01<a></th><td>Statistics about the reads in a BAM.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/FindCorruptedFiles">FindCorruptedFiles<a></th><td>Reads filename from stdin and prints corrupted NGS files (VCF/BAM/FASTQ). </td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFFilterGO">VCFFilterGO<a></th><td>Set the VCF FILTERs on VCF files annotated with SNPEFF or VCP testing wether a Gene belong or not to the descendants of a GO term.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFBed">VCFBed<a></th><td>Annotate a VCF with the content of a BED file indexed with tabix.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar76892">Biostar76892<a></th><td>Fix strand of two paired reads close but on the same strand http://www.biostars.org/p/76892/</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFCompareGT">VCFCompareGT<a></th><td>VCF : compare genotypes of two or more callers for the same samples.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/HowManyBamDict">HowManyBamDict<a></th><td>for @abinouze : quickly find the number of distinct BAM Dictionaries from a set of BAM files.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFPolyX">VCFPolyX<a></th><td>Number of repeated REF bases around POS.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar81455">Biostar81455<a></th><td>Defining precisely the genomic context based on a position http://www.biostars.org/p/81455/</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/NoZeroVariationVCF">NoZeroVariationVCF<a></th><td>creates a VCF containing one fake variation if the input is empty.</td></tr>
</table>


<h3>FixVCF</h3>
<h4>Motivation</h4>
Fix a VCF HEADER when I forget to declare a FILTER or an INFO field in the HEADER
<h4>Compilation</h4>
```bash
ant fixvcf
```
<h4>Options</h4>
<table>
<tr><th>Option</th><th>Description</th></tr>
<tr><td>IN=String</td><td>VCF file/URL to process. Default stdin. </td></tr>
<tr><td>OUT=File</td><td>VCF file to generate. Default stdout.</td></tr>
</table>

<h3>VCFTrio</h3>
<h4>Motivation</h4>
Check for mendelian incompatibilities in a VCF
<h4>Compilation</h4>
```bash
ant vcftrio
```
<h4>Options</h4>
<table>
<tr><th>Option</th><th>Description</th></tr>
<tr><td>PEDIGREE=File</td><td> Pedigree file (plink format) Required. </td></tr>
<tr><td>FILTER=Boolean</td><td>Set filter 'MENDEL' if incompatibilities found.  Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false} </td></tr>
<tr><td>IN=String</td><td>VCF file/URL to process. Default stdin.   Default value: null. </td></tr>
<tr><td>OUT=File</td><td>VCF file to generate. Default stdout.   Default value: null. </td></tr>
</table>

<h3>MapUniProtFeatures</h3>
<h4>Motivation</h4>
map uniprot features on reference genome.

<h4>Options</h4>
<table>
<tr><th>Option</th><th>Description</th></tr>
<tr><td>REF=File</td><td>Reference  Required. </td></tr>
<tr><td>OUT=File</td><td>output name (default: stdout)  Required. </td></tr>
<tr><td>KGURI=String</td><td>KnownGene data  Default value: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz. This option can be set to 'null' to clear the default value. </td></tr>
<tr><td>UNIPROT=String</td><td>Uniprot URL/File  Default value: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz. This option can be set to 'null' to clear the default value. </td></tr>
</table>
<h4>Example</h4>
```bash
$ java  -jar dist/mapuniprot.jar \
	REF=/path/to/human_g1k_v37.fasta \
	UNIPROT=/path/uri/uniprot.org/uniprot_sprot.xml.gz  \
	kgUri=<(curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz" | gunzip -c | awk -F '        ' '{if($2 ~ ".*_.*") next; OFS="       "; gsub(/chr/,"",$2);print;}'   ) |\
	LC_ALL=C sort -t '	' -k1,1 -k2,2n -k3,3n  | uniq | head


1	69090	69144	topological_domain	1000	+	69090	69144	255,0,0	1	54	0
1	69144	69216	transmembrane_region	1000	+	69144	69216	255,0,0	1	72	0
1	69216	69240	topological_domain	1000	+	69216	69240	255,0,0	1	24	0
1	69240	69306	transmembrane_region	1000	+	69240	69306	255,0,0	1	66	0
1	69306	69369	topological_domain	1000	+	69306	69369	255,0,0	1	63	0
1	69357	69636	disulfide_bond	1000	+	69357	69636	255,0,0	1	279	0
1	69369	69429	transmembrane_region	1000	+	69369	69429	255,0,0	1	60	0
1	69429	69486	topological_domain	1000	+	69429	69486	255,0,0	1	57	0
1	69486	69543	transmembrane_region	1000	+	69486	69543	255,0,0	1	57	0
1	69543	69654	topological_domain	1000	+	69543	69654	255,0,0	1	111	0
```
<h3>ExtendBed</h3>
<h4>Motivation</h4>
Extends a BED file by 'X' bases.
<h4>Options</h4>
<table>
<tr><th>Option</th><th>Description</th></tr>
<tr><td>IN=String</td><td>BED Input URI/file. default: stdin  Default value: null. </td></tr>
<tr><td>REF=File</td><td>Reference  Required. </td></tr>
<tr><td>OUT=File</td><td>output name (default: stdout)</td></tr>
<tr><td>EXTEND=Integer</td><td>extend by 'X' bases.  Default value: 0. This option can be set to 'null' to clear the default value. </td></tr>
</table>
<h4>Example</h4>
```bash
head test.bed |\
	java -jar dist/extendbed.jar \
	X=100 REF=human_g1k_v37.fa
```
<h3>Biostar77288</h3>
<h4>Motivation</h4>
Low resolution sequence alignment visualization http://www.biostars.org/p/77288/
<h4>Compilation</h4>
```bash
ant biostar77288
 ```
<h4>Options</h4>
<table>
<tr><th>Option</th><th>Description</th></tr>
<tr><td>IN=String</td><td>MSA file/URI (default:stdin)</td></tr>
<tr><td>ALN_WIDTH=Integer</td><td>Alignment width  Default value: 1000. This option can be set to 'null' to clear the default value. </td></tr>
<tr><td>SEQLOGO=Boolean</td><td>Input is seqLogo (see https://github.com/lindenb/jvarkit#sam4weblogo)  Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false} </td></tr>
</table>
<h4>Example</h4>
```bash
curl -s "http://www.tcoffee.org/Courses/Exercises/saragosa_pb_2010/practicals/practical_2/ex.1.19/file/clustalw.msa" |\
	java -jar dist/biostar77288.jar  > result.svg
```
![ScreenShot](https://raw.github.com/lindenb/jvarkit/master/doc/biostar77288.png)
```bash
$ java -jar dist/sam4weblogo.jar IN=in.bam   REGION="1:630-719" |\
	java -jar dist/biostar77288.jar  SEQLOGO=true > result.svg
```


<h3>VCFBigWig</h3>
<h4>Motivation</h4>
Annotate a VCF with the data of a bigwig file.
<h4>Compilation</h4>
Compiling requires the bigwig library http://code.google.com/p/bigwig/

```bash

$ more build.properties 
(...)
bigwig.dir=/path/to/bigwig
(...)

$ ant vcfbigwig
```
<h4>Options</h4>
<table>
<tr><th>Option</th><th>Description</th></tr>
<tr><td>BW=File</td><td> Path to the bigwig file. The chromosome must have the same names than in the VCF. Required.</td></tr>
<tr><td>INFOTAG=String</td><td>name of the INFO tag. default: name of the bigwig.  Default value: null</td></tr>
<tr><td>IN=String</td><td>VCF file/URL to process. Default stdin.   Default value: null. </td></tr>
<tr><td>OUT=File</td><td>VCF file to generate. Default stdout.   Default value: null. </td></tr>
</table>
<h4>Example</h4>
```bash
 java -jar dist/vcfbigwig.jar \
 	IN=input.vcf.gz \
 	INFOTAG=GERP \
 	BW=gerp.bw
	
##INFO=<ID=GERP,Number=1,Type=Float,Description="Values from bigwig file: com.github.lindenb.jvarkit.tools.vcfbigwig.VCFBigWig BIGWIG=gerp.bw TAG=GERP IN=input.vcf.gz    VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO(...)
A	33926	.	G	A	182	.	GERP=-6.35(...)
A	45365	.	A	G	222	.	GERP=-3.55(...)
```

<h3>VCFTabixml</h3>
<h4>Motivation</h4>
Annotate a value from a vcf+xml file.4th column of the BED indexed with TABIX is a XML string.It will be processed with the xslt-stylesheet and should procuce a xml result &lt;properties>&lt;property key='key1'>value1&lt;/property>&lt;property key='key2'>values1&lt;/property>&lt;/properies> INFO fields. Carriage returns will be removed.Parameters to be passed to the stylesheet: vcfchrom (string) vcfpos(int) vcfref(string) vcfalt(string). Version: 1.0
See also evs2bed

<h4>Compilation</h4>
```bash
ant vcftabixml
```
<h4>Options</h4>
<table>
<tr><th>Option</th><th>Description</th></tr>
<tr><td>BEDFILE=String</td><td> BED file indexed with tabix. The 4th column *is* a XML string.)  Required. </td></tr>
<tr><td>STYLESHEET=File</td><td>x xslt-stylesheet. REQUIRED. Should produce a valid set of INFO field.  Required. </td></tr>
<tr><td>TAGS=File</td><td>file containing extra INFO headers line to add version: 4.1  Required. </td></tr>
<tr><td>IN=String</td><td>VCF file/URL to process. Default stdin.   Default value: null. </td></tr>
<tr><td>OUT=File</td><td>VCF file to generate. Default stdout.   Default value: null. </td></tr>
</table>
<h4>Example</h4>
The header of the tabix-index BED+XML file:
```
1	69427	69428	<snpList><positionString>1:69428</positionString><chrPosition>69428</chrPosition><alleles>G/T</alleles><ua..
1	69475	69476	<snpList><positionString>1:69476</positionString><chrPosition>69476</chrPosition><alleles>C/T</alleles><ua..
1	69495	69496	<snpList><positionString>1:69496</positionString><chrPosition>69496</chrPosition><alleles>A/G</alleles><ua..
1	69510	69511	<snpList><positionString>1:69511</positionString><chrPosition>69511</chrPosition><alleles>G/A</alleles><ua...
1	69589	69590	<snpList><positionString>1:69590</positionString><chrPosition>69590</chrPosition><alleles>A/T</alleles><ua...
1	69593	69594	<snpList><positionString>1:69594</positionString><chrPosition>69594</chrPosition><alleles>C/T</alleles><ua..
1	69619	69620	<snpList><positionString>1:69620</positionString><chrPosition>69620</chrPosition><alleles>T/TA</alleles><u...
1	69744	69745	<snpList><positionString>1:69745</positionString><chrPosition>69745</chrPosition><alleles>CA/C</alleles><u..
1	69760	69761	<snpList><positionString>1:69761</positionString><chrPosition>69761</chrPosition><alleles>T/A</alleles><ua...
1	801942	801943	<snpList><positionString>1:801943</positionString><chrPosition>801943</chrPosition><alleles>T/C</alleles...
(...)
```
The content of the tags.txt file:
```
##INFO=<ID=EVS_AAMAF,Number=1,Type=String,Description="EVS_AAMAF from EVS">
##INFO=<ID=EVS_AVGSAMPLEREADDEPTH,Number=1,Type=String,Description="EVS_AVGSAMPLEREADDEPTH from EVS">
##INFO=<ID=EVS_CLINICALLINK,Number=1,Type=String,Description="EVS_CLINICALLINK from EVS">
##INFO=<ID=EVS_CONFLICTALT,Number=1,Type=String,Description="EVS_CONFLICTALT from EVS">
##INFO=<ID=EVS_CONSERVATIONSCOREGERP,Number=1,Type=String,Description="EVS_CONSERVATIONSCOREGERP from EVS">
##INFO=<ID=EVS_CONSERVATIONSCORE,Number=1,Type=String,Description="EVS_CONSERVATIONSCORE from EVS">
##INFO=<ID=EVS_GENELIST,Number=1,Type=String,Description="EVS_GENELIST from EVS">
##INFO=<ID=EVS_GWASPUBMEDIDS,Number=1,Type=String,Description="EVS_GWASPUBMEDIDS from EVS">
##INFO=<ID=EVS_ONEXOMECHIP,Number=1,Type=String,Description="EVS_ONEXOMECHIP from EVS">
##INFO=<ID=EVS_RSIDS,Number=1,Type=String,Description="EVS_RSIDS from EVS">
##INFO=<ID=EVS_TOTALMAF,Number=1,Type=String,Description="EVS_TOTALMAF from EVS">
##INFO=<ID=EVS_UAMAF,Number=1,Type=String,Description="EVS_UAMAF from EVS">
```
The XSLT file:

```xml
<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0"> 
<xsl:output method="xml"/>
<xsl:param name="vcfchrom"/>
<xsl:param name="vcfpos"/>
<xsl:param name="vcfref"/>
<xsl:param name="vcfalt"/>
<xsl:template match="/">
<properties>
<xsl:apply-templates select="evsData|snpList"/>
</properties>
</xsl:template>
<xsl:template match="evsData">
<xsl:apply-templates select="snpList"/>
</xsl:template>
<xsl:template match="snpList">
<xsl:choose>
 <xsl:when test="chromosome=$vcfchrom and chrPosition=$vcfpos and refAllele=$vcfref">
<xsl:apply-templates select="clinicalLink|rsIds|uaMAF|aaMAF|totalMAF|avgSampleReadDepth|geneList|conservationScore|conservationScoreGERP|gwasPubmedIds|onExomeChip|gwasPubmedIds"/>
  <xsl:if test="altAlleles!=$vcfalt">
  	<entry key="EVS_CONFLICTALT">
  		<xsl:value-of select="altAlleles"/>
  	</entry>		
  	</xsl:if>
  </xsl:when>
  <xsl:otherwise/>
</xsl:choose>
</xsl:template>
<xsl:template match="clinicalLink|rsIds|uaMAF|aaMAF|totalMAF|avgSampleReadDepth|geneList|conservationScore|conservationScoreGERP|onExomeChip|gwasPubmedIds">
<xsl:if test="string-length(normalize-space(text()))&gt;0">
<entry>
<xsl:attribute name="key"><xsl:value-of select="concat('EVS_',translate(name(.),'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ'))"/></xsl:attribute>
<xsl:value-of select="text()"/>
</entry>
</xsl:if>
</xsl:template>
</xsl:stylesheet>
```

Running:

```bash
java -jar dist/vcftabixml.jar \
	IN=input.vcf.gz \
	XSL=evs2vcf.xsl \
	BEDFILE=evs.data.gz \
	TAGS=tags.txt | grep -v "##" | head

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	U0925
A	613326	.	G	A	182	.	AC1=1;AF1=0.5;DP=327;DP4=0,147,0,178;FQ=163;MQ=47;PV4=1,0.21,1.2e-35,1;RPB=6.999638e+00;VDB=8.105711e-05	GT:PL:DP:GQ	0/1:212,0,190:325:99
A	614565	.	A	G	222	.	AC1=2;AF1=1;DP=910;DP4=0,1,71,836;EVS_AAMAF=23.2659;EVS_AVGSAMPLEREADDEPTH=70;EVS_CLINICALLINK=unknown;EVS_CONSERVATIONSCORE=0.0;EVS_CONSERVATIONSCOREGERP=-3.5;EVS_GENELIST=KCNAB2;EVS_GWASPUBMEDIDS=unknown;EVS_ONEXOMECHIP=false;EVS_RSIDS=rs26;EVS_TOTALMAF=30.5519;EVS_UAMAF=33.7209;FQ=-282;MQ=41;PV4=1,0.26,0.032,1;RPB=1.705346e+00;VDB=1.656917e-35	GT:PL:DP:GQ	1/1:255,255,0:908:99
A	614379	.	C	T	225	.	AC1=1;AF1=0.5;DP=979;DP4=33,440,37,456;EVS_AAMAF=2.8179;EVS_AVGSAMPLEREADDEPTH=59;EVS_CLINICALLINK=unknown;EVS_CONSERVATIONSCORE=0.0;EVS_CONSERVATIONSCOREGERP=-4.7;EVS_GENELIST=KCNAB2;EVS_GWASPUBMEDIDS=unknown;EVS_ONEXOMECHIP=false;EVS_RSIDS=rs249;EVS_TOTALMAF=8.8261;EVS_UAMAF=11.4393;FQ=225;MQ=42;PV4=0.8,1,3.8e-152,1;RPB=1.317662e+01;VDB=3.857882e-21	GT:PL:DP:GQ	0/1:255,0,255:966:99
A	614565	.	T	G	222	.	AC1=2;AF1=1;DP=209;DP4=0,0,0,187;FQ=-282;MQ=46;VDB=2.410569e-12	GT:PL:DP:GQ	1/1:255,255,0:187:99
A	614953	.	C	T	225	.	AC1=1;AF1=0.5;DP=810;DP4=197,214,194,195;EVS_AAMAF=2.4603;EVS_AVGSAMPLEREADDEPTH=45;EVS_CLINICALLINK=unknown;EVS_CONSERVATIONSCORE=0.0;EVS_CONSERVATIONSCOREGERP=-3.1;EVS_GENELIST=KCNAB2;EVS_GWASPUBMEDIDS=unknown;EVS_ONEXOMECHIP=false;EVS_RSIDS=rs733;EVS_TOTALMAF=7.0506;EVS_UAMAF=9.4205;FQ=225;MQ=47;PV4=0.62,1,9.2e-58,0.015;RPB=4.560133e+00;VDB=6.880316e-03	GT:PL:DP:GQ	0/1:255,0,255:800:99
A	614922	.	G	A	130	.	AC1=1;AF1=0.5;DP=183;DP4=0,94,0,86;EVS_AAMAF=2.3053;EVS_AVGSAMPLEREADDEPTH=43;EVS_CLINICALLINK=unknown;EVS_CONSERVATIONSCORE=0.0;EVS_CONSERVATIONSCOREGERP=-3.0;EVS_GENELIST=KCNAB2;EVS_GWASPUBMEDIDS=unknown;EVS_ONEXOMECHIP=false;EVS_RSIDS=rs202;EVS_TOTALMAF=7.1215;EVS_UAMAF=9.6386;FQ=133;MQ=44;PV4=1,0.0065,3.7e-51,1;RPB=3.500959e+00;VDB=9.915205e-29	GT:PL:DP:GQ	0/1:160,0,184:180:99
A	614986	.	G	C	188	.	AC1=2;AF1=1;DP=176;DP4=0,0,0,175;FQ=-282;MQ=46;VDB=0.000000e+00	GT:PL:DP:GQ	1/1:221,255,0:175:99
A	615009	.	T	A	125	.	AC1=1;AF1=0.5;DP=103;DP4=45,0,56,0;FQ=120;MQ=45;PV4=1,0.14,2.8e-19,1;RPB=1.520268e+00;VDB=1.539079e-06	GT:PL:DP:GQ	0/1:155,0,148:101:99
A	615037	.	C	T	161	.	AC1=1;AF1=0.5;DP=353;DP4=0,164,0,165;FQ=110;MQ=48;PV4=1,1,1.1e-23,1;RPB=5.549816e+00;VDB=1.486773e-11	GT:PL:DP:GQ	0/1:191,0,137:329:99
```



### FindCorruptedFiles


### NgsFilesSummary

#### Motivation

Scan folders and generate a summary of the files (SAMPLE/BAM SAMPLE/VCF etc..)

#### Compilation

```bash
ant ngsfilessummary
```

#### Options

<table>
<tr><th>Option</th><th>Description</th></tr>
<tr><td>IN=File</td><td>File(s) and/or directories  This option may be specified 0 or more times. </td></tr>
</table>

#### Example

```bash
java -jar dist/ngsfilessummary.jar \
	I=/projects/align01/ \
	VALIDATION_STRINGENCY=SILENT 2> /dev/null

SAMPLE1	BAM	/projects/align01/Samples/SAMPLE1/BAM/SAMPLE1_final.bam	321262321	Wed Jun 26 10:30:07 CEST 2013
SAMPLE1	FASTQ	/projects/data
SAMPLE1	VCF	/projects/align01/Samples/SAMPLE1/VCF/SAMPLE1_variations.freebayes.vcf.gz	184191	Mon Jun 17 14:47:22 CEST 2013
SAMPLE1	VCF	/projects/align01/Samples/SAMPLE1/VCF/SAMPLE1_variations.gatk.vcf.gz	113341	Mon Jun 17 11:57:19 CEST 2013
SAMPLE1	VCF	/projects/align01/Samples/SAMPLE1/VCF/SAMPLE1_variations.samtools.vcf.gz	57518	Mon Jun 17 11:58:49 CEST 2013
SAMPLE2	BAM	/projects/align01/Samples/SAMPLE2/BAM/SAMPLE2_final.bam	286100773	Wed Jun 26 10:47:09 CEST 2013
SAMPLE2	FASTQ	/projects/data
SAMPLE2	VCF	/projects/align01/Samples/SAMPLE2/VCF/SAMPLE2_variations.freebayes.vcf.gz	172970	Mon Jun 17 14:45:51 CEST 2013
SAMPLE2	VCF	/projects/align01/Samples/SAMPLE2/VCF/SAMPLE2_variations.gatk.vcf.gz	106390	Mon Jun 17 11:57:19 CEST 2013
SAMPLE2	VCF	/projects/align01/Samples/SAMPLE2/VCF/SAMPLE2_variations.samtools.vcf.gz	52709	Mon Jun 17 11:58:04 CEST 2013
```

### Biostar77828

#### Motivation

Divide the human genome among X cores, taking into account gaps See http://www.biostars.org/p/77828/

#### Compilation
```bash
ant biostar77828
```
<h4>Options</h4>
<table>
<tr><th>Option</th><th>Description</th></tr>
<tr><td>IN=File</td><td>Bed file input (default stdin)  Default value: null. </td></tr>
<tr><td>MIN_CORE=Integer</td><td>min_core  Default value: 20. This option can be set to 'null' to clear the default value. </td></tr>
<tr><td>MAX_CORE=Integer</td><td>max_core  Default value: 30. This option can be set to 'null' to clear the default value. </td></tr>
<tr><td>N_ITERATIONS=Long</td><td>number of iterations  Default value: 1000000. This option can be set to 'null' to clear the default value. </td></tr>
</table>

### Biostar78285

#### Motivation

Extract regions of genome that have 0 coverage See http://www.biostars.org/p/78285/

#### Compilation
```bash
ant biostar78285
```
<h4>Options</h4>
<table>
<tr><th>Option</th><th>Description</th></tr>
<tr><td>IN=File</td><td>BAM file (sorted on coordinate). Default:stdin  Default value: null. </td></tr>
<tr><td>USECIGAR=Boolean</td><td>scan the CIGAR string & detect the gaps in the reads. Slower & requires more memory  Default value: false.</td></tr>
</table>

#### Example:
```bash
 $ java -jar dist/biostar78285.jar \
 	I=sorted.bam \
 	USECIGAR=false \
 	VALIDATION_STRINGENCY=LENIENT

seq1	1569	1575
seq2	1567	1584
```

### Biostar78400

#### Motivation

add the read group info to the sam file on a per lane basis

#### Compilation
```bash
ant biostar78400
```
<h4>Options</h4>
<table>
<tr><th>Option</th><th>Description</th></tr>
<tr><td>IN=File</td><td>BAM file to process (or stdin).  Default value: null. </td></tr>
<tr><td>OUT=File</td><td>BAM file (or stdout).  Default value: null. </td></tr>
<tr><td>XML=File</td><td>XML desfription of the groups. See below Required. </td></tr>
</table>

#### XML config
the XML should look like this:

```xml
<read-groups>
<flowcell name="HS2000-1259_127">
 <lane index="1">
   <group ID="X1">
     <library>L1</library>
     <platform>P1</platform>
     <sample>S1</sample>
     <platformunit>PU1</platformunit>
     <center>C1</center>
     <description>blabla</description>
   </group>
 </lane>
</flowcell>
<flowcell name="HS2000-1259_128">
 <lane index="2">
   <group ID="x2">
     <library>L2</library>
     <platform>P2</platform>
     <sample>S2</sample>
     <platformunit>PU1</platformunit>
     <center>C1</center>
     <description>blabla</description>
   </group>
 </lane>
</flowcell>
</read-groups>
```

#### Example:

```bash
$ cat input.sam 
@SQ	SN:ref	LN:45
@SQ	SN:ref2	LN:40
HS2000-1259_127:1:1210:15640:52255	163	ref	7	30	8M4I4M1D3M	=	37	39	
TTAGATAAAGAGGATACTG	*	XX:B:S,12561,2,20,112
HS2000-1259_128:2:1210:15640:52255	0	ref	9	30	1S2I6M1P1I1P1I4M2I	*	0	
0	AAAAGATAAGGGATAAA	*

$java -jar dist/biostar78400.jar \
	XML=groups.xml \
	I=input.sam \
 	VALIDATION_STRINGENCY=LENIENT

@HD	VN:1.4	SO:unsorted
@SQ	SN:ref	LN:45
@SQ	SN:ref2	LN:40
@RG	ID:X1	PL:P1	PU:P1	LB:L1	DS:blabla	SM:S1	CN:C1
@RG	ID:x2	PL:P2	PU:P2	LB:L2	DS:blabla	SM:S2	CN:C1
@PG	ID:Biostar78400	PN:Biostar78400	PP:Biostar78400	VN:1.0	(...)
HS2000-1259_127:1:1210:15640:52255	163	ref	7	30	8M4I4M1D3M	=	37	39	TTAGATAAAGAGGATACTG	*	RG:Z:X1	XX:B:S,12561,2,20,112
HS2000-1259_128:2:1210:15640:52255	0	ref	9	30	1S2I6M1P1I1P1I4M2I	*	0	0AAAAGATAAGGGATAAA	*	RG:Z:x2
```

### BamStats04

#### Motivation

 Coverage statistics for a BED file. It uses the Cigar string instead of the start/end to compute the coverage
 
#### Compilation
```bash
ant bamstats04
```
<h4>Options</h4>
<table>
<tr><th>Option</th><th>Description</th></tr>
<tr><td>IN=File</td><td>BAM file to process.  Required. </td></tr>
<tr><td>BEDILE=File</td><td>BED File.  Required. </td></tr>
<tr><td>NO_DUP=Boolean</td><td>discard duplicates  Default value: true. This option can be set to 'null' to clear the default value. Possible values: {true, false} </td></tr>
<tr><td>NO_ORPHAN=Boolean</td><td>discard not properly paired  Default value: true. This option can be set to 'null' to clear the default value. Possible values: {true, false} </td></tr>
<tr><td>NO_VENDOR=Boolean</td><td>discard failing Vendor Quality  Default value: true. This option can be set to 'null' to clear the default value. Possible values: {true, false} </td></tr>
<tr><td>MMQ=Integer</td><td>min mapping quality  Default value: 0. This option can be set to 'null' to clear the default value. </td></tr>
<tr><td>MIN_COVERAGE=Integer</td><td>min coverage to say the position is not covered  Default value: 0. This option can be set to 'null' to clear the default value. </td></tr>
</table>


#### Example:

```bash
$ java -jar dist/bamstats04.jar \
	BED=data.bed \
	I=f.bam
#chrom	start	end	length	mincov	maxcov	mean	nocoveragepb	percentcovered
1	429665	429785	120	42	105	72.36666666666666	0	100
1	430108	430144	36	9	9	9.0	0	100
1	439811	439904	93	0	36	3.6451612903225805	21	77
1	550198	550246	48	1325	1358	1344.4583333333333	0	100
1	629855	629906	51	223	520	420.70588235294116	0	100
1	689960	690029	69	926	1413	1248.9420289855072	0	100
1	690852	690972	120	126	193	171.24166666666667	0	100
1	787283	787406	123	212	489	333.9756097560976	0	100
1	789740	789877	137	245	688	528.6715328467153	0	1
```

### VCFAnnoBam

#### Motivation

Annotate a VCF with the Coverage statistics of a BAM file +  BED file of capture.
It uses the Cigar string instead of the start/end to get the voverage

 
#### Compilation  #### 
```bash
ant vcfannobam
```
<h4>Options</h4>
<table>
<tr><th>Option</th><th>Description</th></tr>
<tr><td>BEDILE=File</td><td>BED File capture.  Required. </td></tr>
<tr><td>BAMFILE=File</td><td>indexed BAM File.  This option must be specified at least 1 times. </td></tr>
<tr><td>MMQ=Integer</td><td>min mapping quality  Default value: 0. This option can be set to 'null' to clear the default value. </td></tr>
<tr><td>MIN_COVERAGE=Integer</td><td>min coverage to say the position is not covered  Default value: 0. This option can be set to 'null' to clear the default value. </td></tr>
<tr><td>IN=String</td><td>VCF file/URL to process. Default stdin.   Default value: null. </td></tr>
<tr><td>OUT=File</td><td>VCF file to generate. Default stdout.   Default value: null. </td></tr>
</table>


Example:  #### 

```bash
$  java -jar dist/vcfannobam.jar \
		BAM=input.bam\
		IN=input.vcf.gz \
		BED=capture.bed \
		VALIDATION_STRINGENCY=SILENT

(...)
##INFO=<ID=CAPTURE,Number=1,Type=String,Description="Capture stats: Format is (start|end|mean|min|max|length|not_covered|percent_covered) ">
(...)
2	16100665	.	A	T	13261.77	.	CAPTURE=16100619|16100715|1331.96|1026.0|1773.0|97|0|100
2	178395141	.	T	A	1940.77	.	CAPTURE=178394991|178395199|193.11|100.0|276.0|209|0|100
(...)
```


### EVS2Bed ###

#### Motivation ####

Download data from EVS http://evs.gs.washington.edu/EVS as a BED chrom/start/end/XML
For later use, see vcftabixml

#### Compilation  #### 

```bash
ant evs2bed
```	

#### Example  ####
```bash
$ java  -jar dist/evs2bed.jar L=10 2> /dev/null | cut -c 1-100
1	69427	69428	<snpList><positionString>1:69428</positionString><chrPosition>69428</chrPosition><alle..
1	69475	69476	<snpList><positionString>1:69476</positionString><chrPosition>69476</chrPosition><alle..
1	69495	69496	<snpList><positionString>1:69496</positionString><chrPosition>69496</chrPosition><alle..
1	69510	69511	<snpList><positionString>1:69511</positionString><chrPosition>69511</chrPosition><alle..
```
