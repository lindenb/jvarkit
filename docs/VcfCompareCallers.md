# VcfCompareCallers


## Usage

```
Usage: vcfcomparecallers [options] Files
  Options:
    -B, --bed
      Limit to variants in that BED region
    -e, --examplefile
      Write a few Variants in this XML file. Optional
    -h, --help
      print help and exits
    -c, --homref2nocall
      Treat HomRef as No Call (created when comparing merged vcf with GATK: 
      there is no homref, everything is nocall)
      Default: false
    -d, --noindel
      Ignore indels (SNP only)
      Default: false
    -n, --num
      number of variants to dump in the example file
      Default: 10
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exits
    --vertical
      Vertical layout, to use with group-by tools
      Default: false

```


## Description

Compare two VCFs and print common/exclusive information for each sample/genotype

## Compilation

### Requirements / Dependencies

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make vcfcomparecallers
```

The *.jar libraries are not included in the main jar file, so you shouldn't move them (https://github.com/lindenb/jvarkit/issues/15#issuecomment-140099011 ).
The required libraries will be downloaded and installed in the `dist` directory.

### edit 'local.mk' (optional)

The a file **local.mk** can be created edited to override/add some definitions.

For example it can be used to set the HTTP proxy:

```
http.proxy.host=your.host.com
http.proxy.port=124567
```
## Source code 

https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcmp/VcfCompareCallers.java

## Contribute

- Issue Tracker: http://github.com/lindenb/jvarkit/issues
- Source Code: http://github.com/lindenb/jvarkit

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfcomparecallers** ? https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md

The current reference is:

http://dx.doi.org/10.6084/m9.figshare.1425030

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> http://dx.doi.org/10.6084/m9.figshare.1425030






### Synopsis




```
$ java -jar dist/vcfcomparecallers.jar file1.vcf(.gz) stdin 
$ java -jar dist/vcfcomparecallers.jar file1.vcf(.gz) file2.vcf(.gz) 

```



both vcf **must** be sorted on CHROM/POS/REF



### Column headers


 *  off_target_only_1 : variant for this sample is found in vcf1 but not in vcf2, but variant is off target
 *  off_target_only_2 : variant for this sample is found in vcf2 but not in vcf1, but variant is off target
 *  off_target_both : variant for this sample is found in both vcf, but it variant off target
 *  unique_to_file_1 :  variant for this sample is found in vcf1 but not in vcf2 (in target)
 *  unique_to_file_1_snp : variant for this sample is found in vcf1 but not in vcf2 (in target) and it is a snp
 *  unique_to_file_1_indel : variant for this sample is found in vcf1 but not in vcf2 (in target) and it is an indel
 *  unique_to_file_2 :  variant for this sample is found in vcf2 but not in vcf1 (in target)
 *  unique_to_file_2_snp: variant for this sample is found in vcf2 but not in vcf1 (in target) and it is a snp
 *  unique_to_file_2_indel: variant for this sample is found in vcf2 but not in vcf1 (in target) and it is an indel
 *  both_missing: there is no genotype for this variant (chrom/pos/ref), while some other samples can have a called genotype.
 *  common_context: a genotype is available for this sample for this context in vcf1 and vcf2.
 *  common_context_snp:  a genotype is available for this sample for this context in vcf1 and vcf2 and it's a snp.
 *  common_context_indel:  a genotype is available for this sample for this context in vcf1 and vcf2 and it's an indel.
 *  common_context_discordant_id:  a genotype is available for this sample for this context in vcf1 and vcf2  but the column ID was not the same.
 *  common_context_discordant_filter:  a genotype is available for this sample for this context in vcf1 and vcf2  but the column FILTER is set in one and not in the other
 *  called_and_same : vcf1 and vcf2 have the same genotype for this sample and variant (chrom/pos/ref).
 *  called_and_same_hom_ref :  vcf1 and vcf2 have the same hom-ref genotype for this sample and variant (chrom/pos/ref).
 *  called_and_same_hom_var :  vcf1 and vcf2 have the same hom-var genotype for this sample and variant (chrom/pos/ref).
 *  called_and_same_het : vcf1 and vcf2 have the same het genotype for this sample and variant (chrom/pos/ref).
 *  called_but_discordant : vcf1 and vcf2 don't have the same genotype for this variant and sample.
 *  called_but_discordant_snp : vcf1 and vcf2 don't have the same genotype for this variant and sample and the variant is a SNP.
 *  called_but_discordant_hom1_het2 :  vcf1 and vcf2 don't have the same genotype for this variant and sample. vcf1 is hom and vcf2 is het
 *  called_but_discordant_het1_hom2 :  vcf1 and vcf2 don't have the same genotype for this variant and sample. vcf1 is het and vcf2 is hom
 *  called_but_discordant_hom1_hom2 :  vcf1 and vcf2 don't have the same genotype for this variant and sample. vcf1 is hom and vcf2 is hom
 *  called_but_discordant_het1_het2 :  vcf1 and vcf2 don't have the same genotype for this variant and sample. vcf1 is het and vcf2 is het
 *  called_but_discordant_others :  vcf1 and vcf2 don't have the same genotype for this variant and sample. other cases.







### Example



```
$ java -jar dist-1.128/vcfcomparecallers.jar  Proj1.samtools.vcf.gz  Proj1.varscan.vcf.gz
#Sample	unique_to_file_1	unique_to_file_1_snp	unique_to_file_1_indel	unique_to_file_2	unique_to_file_2_snp	unique_to_file_2_indel	both_missing	common_context	common_context_snp	common_context_indel	common_context_discordant_id	called_and_same	called_and_same_hom_ref	called_and_same_hom_varcalled_and_same_het	called_but_discordant	called_but_discordant_hom1_het2called_but_discordant_het1_hom2	called_but_discordant_hom1_hom2	called_but_discordant_het1_het2	called_but_discordant_others
B00G5XG	43739	15531	27518	0	10773	11730	2182	558753	535010	22508	55052	1043356	0	26920	41136	3047	698	1993	152	204	0
B00G74M	43629	15445	27503	0	10739	11747	2346	558716	534939	22526	55092	1043355	0	27962	40295	2910	742	1823	164	181	0
B00G5XF	43542	15344	27515	0	10742	11691	2185	559017	535236	22533	55089	1044311	0	26842	40961	2960	809	1821	157	173	0
B00G74L	43705	15461	27543	0	10765	11745	2356	558606	534872	22509	55053	1041955	0	26849	42430	2989	725	1904	175	185	0
B00G5XE	43589	15393	27515	0	10764	11708	2425	558691	534970	22481	55052	1042648	0	27088	41698	2974	746	1906	152	170	0

```





```

$ java -jar dist-1.128/vcfcomparecallers.jar  Proj1.samtools.vcf.gz  Proj1.varscan.vcf.gz | verticalize

>>> 2
$1	#Sample	B00G5XG
$2	unique_to_file_1	43739
$3	unique_to_file_1_snp	15531
$4	unique_to_file_1_indel	27518
$5	unique_to_file_2	0
$6	unique_to_file_2_snp	10773
$7	unique_to_file_2_indel	11730
$8	both_missing	2182
$9	common_context	558753
$10	common_context_snp	535010
$11	common_context_indel	22508
$12	common_context_discordant_id	55052
$13	called_and_same	1043356
$14	called_and_same_hom_ref	0
$15	called_and_same_hom_var	26920
$16	called_and_same_het	41136
$17	called_but_discordant	3047
$18	called_but_discordant_hom1_het2	698
$19	called_but_discordant_het1_hom2	1993
$20	called_but_discordant_hom1_hom2	152
$21	called_but_discordant_het1_het2	204
$22	called_but_discordant_others	0
<<< 2

>>> 3
$1	#Sample	B00G74M
$2	unique_to_file_1	43629
$3	unique_to_file_1_snp	15445
$4	unique_to_file_1_indel	27503
$5	unique_to_file_2	0
$6	unique_to_file_2_snp	10739
$7	unique_to_file_2_indel	11747
$8	both_missing	2346
$9	common_context	558716
$10	common_context_snp	534939
$11	common_context_indel	22526
$12	common_context_discordant_id	55092
$13	called_and_same	1043355
$14	called_and_same_hom_ref	0
$15	called_and_same_hom_var	27962
$16	called_and_same_het	40295
$17	called_but_discordant	2910
$18	called_but_discordant_hom1_het2	742
$19	called_but_discordant_het1_hom2	1823
$20	called_but_discordant_hom1_hom2	164
$21	called_but_discordant_het1_het2	181
$22	called_but_discordant_others	0
<<< 3

>>> 4
$1	#Sample	B00G5XF
$2	unique_to_file_1	43542
$3	unique_to_file_1_snp	15344
$4	unique_to_file_1_indel	27515
$5	unique_to_file_2	0
$6	unique_to_file_2_snp	10742
$7	unique_to_file_2_indel	11691
$8	both_missing	2185
$9	common_context	559017
$10	common_context_snp	535236
$11	common_context_indel	22533
$12	common_context_discordant_id	55089
$13	called_and_same	1044311
$14	called_and_same_hom_ref	0
$15	called_and_same_hom_var	26842
$16	called_and_same_het	40961
$17	called_but_discordant	2960
$18	called_but_discordant_hom1_het2	809
$19	called_but_discordant_het1_hom2	1821
$20	called_but_discordant_hom1_hom2	157
$21	called_but_discordant_het1_het2	173
$22	called_but_discordant_others	0
<<< 4


<<< 4


```





### Example


the following XSLT stylesheet can be used to produce a HTML table for a few differences:



```

<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns="http://www.w3.org/1999/xhtml" version="1.0">
  <xsl:output method="xml"/>
  <xsl:template match="/">
  	<table>
  	<thead>
  		<tr><th>Category</th><th>Sample</th><th>Variant 1</th><th>Genotype 1</th><th>Variant 2</th><th>Genotype 2</th></tr>
  	</thead>
  	<tbody>
    <xsl:apply-templates select="compare-callers/diff">
    	 <xsl:sort select="@sample" />
    	 <xsl:sort select="@type" /> 
    </xsl:apply-templates>
    </tbody>
    </table>
  </xsl:template>
  
  <xsl:template match="diff">
  	<tr>
  		<td><xsl:value-of select="@type"/></td>
  		<td><xsl:value-of select="@sample"/></td>
  		<td><xsl:apply-templates select="variant[@file='1']"/></td>
  		<td><xsl:apply-templates select="variant[@file='1']/genotype"/></td>
  		<td><xsl:apply-templates select="variant[@file='2']"/></td>
  		<td><xsl:apply-templates select="variant[@file='2']/genotype"/></td>
  	</tr>
  	<xsl:text>
</xsl:text>
  </xsl:template>
 
   <xsl:template match="variant">
    <xsl:value-of select="concat('(',@type,') ',chrom,':',pos,' ',id,' ',ref,'/',alts)"/>
   </xsl:template>
 
  <xsl:template match="genotype">
  	<xsl:value-of select="concat('(',@type,')')"/>
  	<xsl:for-each select="allele">
  	<xsl:if test="position()>1">/</xsl:if>
  	<xsl:value-of select="."/>
  	</xsl:for-each>
  	<xsl:if test="dp"> DP:<xsl:value-of select="dp"/></xsl:if>
   </xsl:template>
 

 
</xsl:stylesheet>

```



This other stylesheet can be used to print the context of the variant for a given BAM (pipe the output to bash )



```

<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns="http://www.w3.org/1999/xhtml" version="1.0">

  <xsl:output method="text"/>
  <xsl:template match="/">

    <xsl:apply-templates select="compare-callers/diff">
    	 <xsl:sort select="@sample" />
    	 <xsl:sort select="@type" /> 
    </xsl:apply-templates>
  </xsl:template>
  
  <xsl:template match="diff">

  	
  	
  	<xsl:variable name="type" select="@type"/>  	
  	<xsl:variable name="sample" select="@sample"/>  	

  	
  	echo "Sample **<xsl:value-of select="@sample"/>** and category **<xsl:value-of select="@type"/>** at **<xsl:choose>
  		<xsl:when test="variant[@file='1']">
  			<xsl:apply-templates select="variant[@file='1']"/>
  		</xsl:when>
  		<xsl:when test="variant[@file='2']">
  			<xsl:apply-templates select="variant[@file='2']"/>
  		</xsl:when>
  	</xsl:choose>**"
  	
  	echo '```'
  	
  	COLUMNS=30 samtools tview -d T -p "<xsl:choose>
  		<xsl:when test="variant[@file='1']">
  			<xsl:apply-templates select="variant[@file='1']" mode="upstream"/>
  		</xsl:when>
  		<xsl:when test="variant[@file='2']">
  			<xsl:apply-templates select="variant[@file='2']" mode="upstream"/>
  		</xsl:when>
  	</xsl:choose>" path/to/file.bam ref.fasta
  	echo '```'
  	

  </xsl:template>
 
   <xsl:template match="variant">
    <xsl:value-of select="concat(chrom,':',pos)"/>
   </xsl:template>
 
    <xsl:template match="variant" mode="upstream">
    <xsl:value-of select="concat(chrom,':',number(pos)-10)"/>
   </xsl:template>
 
</xsl:stylesheet>

```


