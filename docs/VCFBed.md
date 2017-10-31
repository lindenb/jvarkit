# VCFBed

Transfer information from a BED to a VCF


## Usage

```
Usage: vcfbed [options] Files
  Options:
    -B, --bed
      Tribble or Tabix bed file
    -fn, --filternooverlap
      if defined, set this as a FILTER column if not any BED line overlap a 
      variant 
    -fo, --filteroverlap
      if defined, set this as a FILTER column if one or more BED line overlap 
      a variant
    -f, --format
      format pattern ${xx} will be replaced by column xx in the bed line. 
      Empty lines will be ignored (no tag) but the FILTERs will be set.
      Default: ${1}:${2}-${3}
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -ignoreFiltered, --ignoreFiltered
      [20171031] Ignore FILTERed Variants (should be faster)
      Default: false
    -m, --map
      unindexed bed file, will be loaded in memory (faster than tribble/tabix 
      but memory consumming)
    -o, --output
      Output file. Optional . Default: stdout
    -T, --tag
      use the following INFO tag name
      Default: VCFBED
    --version
      print version and exit

```


## Keywords

 * bed
 * vcf
 * annotation



## See also in Biostars

 * [https://www.biostars.org/p/247224](https://www.biostars.org/p/247224)


## Compilation

### Requirements / Dependencies

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) and avoid OpenJdk, use the java from Oracle. Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make vcfbed
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfbed/VCFBed.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfbed/VCFBed.java)


<details>
<summary>Git History</summary>

```
Mon Aug 7 09:53:19 2017 +0200 ; fixed unicode problems after https://github.com/lindenb/jvarkit/issues/82 ; https://github.com/lindenb/jvarkit/commit/68254c69b027a9ce81d8b211447f1c0bf02dc626
Tue Jun 6 18:06:17 2017 +0200 ; postponed vcf ; https://github.com/lindenb/jvarkit/commit/bcd52318caf3cd76ce8662485ffaacaabde97caf
Mon May 15 17:17:02 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/fc77d9c9088e4bc4c0033948eafb0d8e592f13fe
Tue May 9 10:40:20 2017 +0200 ; moving to jcommander ; https://github.com/lindenb/jvarkit/commit/88cfdecb60c1f193ae8b3176ad86181c4a15256b
Fri Apr 14 10:22:16 2017 +0200 ; knime helper ; https://github.com/lindenb/jvarkit/commit/51679edcfb691b8851c06881599c6f1c7a65af34
Thu Apr 13 09:44:05 2017 +0200 ; vcfbed ; https://github.com/lindenb/jvarkit/commit/3e4230782d626f544312689799b48fbf3faf3745
Tue Aug 2 17:12:59 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/f9842fcf255e7b7f3468ca5e2abf947f9c8dd9e3
Wed Mar 23 11:45:49 2016 +0100 ; vcfbed cont ; https://github.com/lindenb/jvarkit/commit/8e7370663fa5f3824ebab9de304a2ee4d90619c1
Sun Feb 28 14:54:44 2016 +0100 ; burden week-end ; https://github.com/lindenb/jvarkit/commit/2f49ec436743934639d0adf51b55a577db7ee3d2
Fri Jun 5 12:42:21 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/cc909f9f4ceea181bb65e4203e3fdbde176c6f2f
Mon Jun 1 15:27:11 2015 +0200 ; change getChrom() to getContig() ; https://github.com/lindenb/jvarkit/commit/5abd60afcdc2d5160164ae6e18087abf66d8fcfe
Fri May 15 18:07:33 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/346bd7d374bfe8f2c969de98ed176060b234f0e1
Thu Apr 30 16:54:43 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/eb6a52869eb3b9b8bf048a079e1d7a96bccab2bf
Wed Apr 29 09:42:54 2015 +0200 ; moved tools vcfpolyx & vcfbed to std argc/argv #tweet ; https://github.com/lindenb/jvarkit/commit/84c8fcb3aa33e74c0a92c6986629a761330c9afe
Tue Apr 28 15:22:51 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/10924ba7a0d78e208109157b92432d6171247c21
Mon May 12 14:06:30 2014 +0200 ; continue moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/011f098b6402da9e204026ee33f3f89d5e0e0355
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Fri Oct 11 15:39:02 2013 +0200 ; picard v.100: deletion of VcfIterator :-( ; https://github.com/lindenb/jvarkit/commit/e88fab449b04aed40c2ff7f9d0cf8c8b6ab14a31
Wed Sep 18 17:39:34 2013 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/7713393b28bc5004c6d8047a24bb9b84d2294fef
Fri Sep 6 15:11:11 2013 +0200 ; moved code for latest version of picard (1.97). Using VCFIterator instead of ASciiLineReader ; https://github.com/lindenb/jvarkit/commit/810877c10406a017fd5a31dacff7e8401089d429
Sun Jul 21 14:17:59 2013 +0200 ; vcf trios, added git HASH in METAINF/Manifest ; https://github.com/lindenb/jvarkit/commit/1854d3695563b91471861164f5e8903042493470
Thu Jul 18 10:23:15 2013 +0200 ; vcf bed ; https://github.com/lindenb/jvarkit/commit/483c2c8b340c33753e6b45dfdc1577803c01135e
Wed Jul 17 17:19:54 2013 +0200 ; many things, updating code, moving sources... ; https://github.com/lindenb/jvarkit/commit/e81829c941f1d0a4ae8334f463ab5cf8ea740c07
Wed Jul 17 08:39:12 2013 +0200 ; tmp ; https://github.com/lindenb/jvarkit/commit/8dc0f2623a4e6e5a5c849ca8083efe689819fea9
Tue Jul 16 19:56:40 2013 +0200 ; disease ontology ; https://github.com/lindenb/jvarkit/commit/3ca2d21908b547d021a7a58e56693956f395a167
Mon Jul 15 12:55:01 2013 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/982384a4b0b74e7c7574f3c897447d7816f0f458
Mon May 6 18:56:46 2013 +0200 ; moving to git ; https://github.com/lindenb/jvarkit/commit/55158d13f0950f16c4a3cc3edb92a87905346ee1
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfbed** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

Map the NCBI biosystems to a BED file using the following script:     https://gist.github.com/6024788 

```
$ gunzip -c ~/ncbibiosystem.bed.gz | head
1	69091	70008	79501	106356	30	Signaling_by_GPCR
1	69091	70008	79501	106383	50	Olfactory_Signaling_Pathway
1	69091	70008	79501	119548	40	GPCR_downstream_signaling
1	69091	70008	79501	477114	30	Signal_Transduction
1	69091	70008	79501	498	40	Olfactory_transduction
1	69091	70008	79501	83087	60	Olfactory_transduction
1	367640	368634	26683	106356	30	Signaling_by_GPCR
1	367640	368634	26683	106383	50	Olfactory_Signaling_Pathway
1	367640	368634	26683	119548	40	GPCR_downstream_signaling
1	367640	368634	26683	477114	30	Signal_Transduction
```

Now, annotate a remote VCF with the data of NCBI biosystems.

```
curl "https://raw.github.com/arq5x/gemini/master/test/test1.snpeff.vcf" |\
 sed 's/^chr//' |\
 java -jar  dist/vcfbed.jar -B ~/ncbibiosystem.bed.gz -T NCBIBIOSYS  -f '($4|$5|$6|$7)' |\
 grep -E '(^#CHR|NCBI)'

##INFO=<ID=NCBIBIOSYS,Number=.,Type=String,Description="metadata added from /home/lindenb/ncbibiosystem.bed.gz . Format was ($4|$5|$6|$7)">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	1094PC0005	1094PC0009	1094PC0012	1094PC0013
1	69270	.	A	G	2694.18	.	AC=40;AF=1.000;AN=40;DP=83;Dels=0.00;EFF=SYNONYMOUS_CODING(LOW|SILENT|tcA/tcG|S60|305|OR4F5|protein_coding|CODING|ENST00000335137|exon_1_69091_70008);FS=0.000;HRun=0;HaplotypeScore=0.0000;InbreedingCoeff=-0.0598;MQ=31.06;MQ0=0;NCBIBIOSYS=(79501|119548|40|GPCR_downstream_signaling),(79501|106356|30|Signaling_by_GPCR),(79501|498|40|Olfactory_transduction),(79501|83087|60|Olfactory_transduction),(79501|477114|30|Signal_Transduction),(79501|106383|50|Olfactory_Signaling_Pathway);QD=32.86	GT:AD:DP:GQ:PL	./.	./.	1/1:0,3:3:9.03:106,9,0	1/1:0,6:6:18.05:203,18,0
1	69511	.	A	G	77777.27	.	AC=49;AF=0.875;AN=56;BaseQRankSum=0.150;DP=2816;DS;Dels=0.00;EFF=NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Aca/Gca|T141A|305|OR4F5|protein_coding|CODING|ENST00000335137|exon_1_69091_70008);FS=21.286;HRun=0;HaplotypeScore=3.8956;InbreedingCoeff=0.0604;MQ=32.32;MQ0=0;MQRankSum=1.653;NCBIBIOSYS=(79501|119548|40|GPCR_downstream_signaling),(79501|106356|30|Signaling_by_GPCR),(79501|498|40|Olfactory_transduction),(79501|83087|60|Olfactory_transduction),(79501|477114|30|Signal_Transduction),(79501|106383|50|Olfactory_Signaling_Pathway);QD=27.68;ReadPosRankSum=2.261	GT:AD:DP:GQ:PL	./.	./.	0/1:2,4:6:15.70:16,0,40	0/1:2,2:4:21.59:22,0,40</h:pre>
```

Another example:

```
$ tabix -h dbsnp138_00-All.vcf.gz "19:58864565-58865165" | sed '/^[^#]/s/^/chr/' |\
java -jar dist/vcfbed.jar -m your.bed -f '${1}|${2}|${3}|${4}&${5}'

##INFO=<ID=VCFBED,Number=.,Type=String,Description="metadata added from your.bed . Format was ${1}|${2}|${3}|${4}&${5}">
(...)
chr19   58864911    rs113760967 T   C   .   .   GNO;OTHERKG;R5;RS=113760967;RSPOS=58864911;SAO=0;SLO;SSR=0;VC=SNV;VCFBED=chr19|58864565|58865165|A1BG&58864865;VP=0x050100020001000102000100;WGT=1;dbSNPBuildID=132
chr19   58865054    rs893183    T   C   .   .   CAF=[0.1299,0.8701];COMMON=1;G5;GNO;HD;KGPROD;KGPhase1;KGPilot123;OTHERKG;PH3;R5;RS=893183;RSPOS=58865054;RV;SAO=0;SLO;SSR=0;VC=SNV;VCFBED=chr19|58864565|58865165|A1BG&58864865;VLD;VP=0x05010002000115051f000100;WGT=1;dbSNPBuildID=86
chr19   58865068    rs893182    T   C   .   .   CAF=[0.1299,0.8701];COMMON=1;G5;GNO;HD;KGPROD;KGPhase1;KGPilot123;OTHERKG;PH3;R5;RS=893182;RSPOS=58865068;RV;SAO=0;SLO;SSR=0;VC=SNV;VCFBED=chr19|58864565|58865165|A1BG&58864865;VLD;VP=0x05010002000115051f000100;WGT=1;dbSNPBuildID=86
chr19   58865082    rs893181    A   T   .   .   CAF=[0.1295,0.8705];COMMON=1;G5;GNO;HD;KGPROD;KGPhase1;KGPilot123;OTHERKG;PH3;R5;RS=893181;RSPOS=58865082;RV;SAO=0;SLO;SSR=0;VC=SNV;VCFBED=chr19|58864565|58865165|A1BG&58864865;VLD;VP=0x05010002000115051f000100;WGT=1;dbSNPBuildID=86
chr19   58865091    rs893180    A   G   .   .   CAF=[0.1299,0.8701];COMMON=1;G5;GNO;HD;KGPROD;KGPhase1;KGPilot123;OTHERKG;R5;RS=893180;RSPOS=58865091;RV;SAO=0;SLO;SSR=0;VC=SNV;VCFBED=chr19|58864565|58865165|A1BG&58864865;VLD;VP=0x05010002000115051e000100;WGT=1;dbSNPBuildID=86
chr19   58865112    rs188818621 C   T   .   .   CAF=[0.9954,0.004591];COMMON=1;KGPROD;KGPhase1;R5;RS=188818621;RSPOS=58865112;SAO=0;SSR=0;VC=SNV;VCFBED=chr19|58864565|58865165|A1BG&58864865;VP=0x050000020001000014000100;WGT=1;dbSNPBuildID=135
chr19   58865164    rs80109863  C   T   .   .   CAF=[0.9949,0.005051];COMMON=1;GNO;KGPROD;KGPhase1;OTHERKG;R5;RS=80109863;RSPOS=58865164;SAO=0;SSR=0;VC=SNV;VCFBED=chr19|58864565|58865165|A1BG&58864865;VP=0x050000020001000116000100;WGT=1;dbSNPBuildID=132
```


