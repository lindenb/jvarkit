# VcfRebase

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Restriction sites overlaping variations in a vcf


## Usage

```
Usage: vcfrefbase [options] Files
  Options:
    -A, --attribute
      VCF INFO attribute
      Default: ENZ
    -E, -enzyme, --enzyme
      restrict to that enzyme name. Default: use all enzymes
      Default: []
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    --outputbcf
      Output bcf (for streams)
      Default: false
    -R, -reference, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --vcfcreateindex
      VCF, create tribble or tabix Index when writing a VCF/BCF to a file.
      Default: false
    --vcfmd5
      VCF, create MD5 checksum when writing a VCF/BCF to a file.
      Default: false
    --version
      print version and exit
    -w, -weight, --weight
      min enzyme weight 6 = 6 cutter like GAATTC, 2 = 2 cutter like ATNNNNNNAT
      Default: 5.0

```


## Keywords

 * vcf
 * rebase
 * restriction
 * enzyme


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfrefbase
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfrebase/VcfRebase.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfrebase/VcfRebase.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfrebase/VcfRebaseTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfrebase/VcfRebaseTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfrefbase** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
 ## Example
 

 ```
 $  curl -s  "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz" |\
   gunzip -c  |\
   java -jar dist/vcfrebase.jar -R  human_g1k_v37.fasta -w 8   |\
   grep -E '(#|ENZ=)' 

##fileformat=VCFv4.1
##INFO=<ID=ENZ,Number=.,Type=String,Description="Enzyme overlapping: Format: (Name,Site,Sequence,pos-1,strand)">
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	256138	rs372711491	T	G	.	.	ENZ=(SwaI|ATTT^AAAT|ATTTAAAT|256139|+),(SwaI|ATTT^AAAT|ATTTAAAt|256131|+);OTHERKG;RS=372711491;RSPOS=256138;SAO=0;SSR=0;VC=SNV;VP=0x050000000001000002000100;WGT=1;dbSNPBuildID=138
1	744076	rs142181107	A	C	.	.	CAF=[0.9775,0.0225];COMMON=1;ENZ=(SwaI|ATTT^AAAT|ATTTAAaT|744070|+);KGPROD;KGPhase1;RS=142181107;RSPOS=744076;SAO=0;SSR=0;VC=SNV;VP=0x050000000001100014000100;WGT=1;dbSNPBuildID=134
1	762592	rs71507462	C	G	.	.	CAF=[0.2755,0.7245];COMMON=1;ENZ=(MauBI|CG^CGCGCG|cGCGCGCG|762592|+);GNO;KGPROD;KGPhase1;OTHERKG;R5;RS=71507462;RSPOS=762592;SAO=0;SLO;SSR=0;VC=SNV;VP=0x050100020001100116000100;WGT=1;dbSNPBuildID=130
1	780347	rs202219272	TTTAA	T	.	.	CAF=[0.9867,0.01331];COMMON=1;ENZ=(PacI|TTAAT^TAA|ttaaTTAA|780348|+);INT;KGPROD;KGPhase1;KGPilot123;OTHERKG;RS=202219272;RSPOS=780348;SAO=0;SSR=0;VC=DIV;VP=0x05000008000110001e000200;WGT=1;dbSNPBuildID=137
1	780361	rs111333032	A	T	.	.	ENZ=(SwaI|ATTT^AAAT|ATTTAAaT|780355|+);GNO;INT;OTHERKG;RS=111333032;RSPOS=780361;SAO=0;SLO;SSR=0;VC=SNV;VP=0x050100080001000102000100;WGT=1;dbSNPBuildID=132
1	786891	rs183914415	C	T	.	.	CAF=[0.9752,0.02479];COMMON=1;ENZ=(AscI|GG^CGCGCC|GGCGCGCC|786892|+);INT;KGPROD;KGPhase1;RS=183914415;RSPOS=786891;SAO=0;SSR=0;VC=SNV;VP=0x050000080001100014000100;WGT=1;dbSNPBuildID=135
1	820332	rs201213736	A	G,T	.	.	ENZ=(SwaI|ATTT^AAAT|ATTTAAaT|820326|+);OTHERKG;RS=201213736;RSPOS=820332;SAO=0;SSR=0;VC=SNV;VP=0x050000000001000002000100;WGT=1;dbSNPBuildID=137
1	820333	rs201439577	T	C	.	.	ENZ=(SwaI|ATTT^AAAT|ATTTAAAt|820326|+);OTHERKG;RS=201439577;RSPOS=820333;SAO=0;SSR=0;VC=SNV;VP=0x050000000001000002000100;WGT=1;dbSNPBuildID=137
1	822560	rs200342299	C	T	.	.	ENZ=(SfiI|GGCCNNNN^NGGCC|GGCCAGcTTGGCC|822554|+);OTHERKG;RS=200342299;RSPOS=822560;SAO=0;SSR=0;VC=SNV;VP=0x050000000001000002000100;WGT=1;dbSNPBuildID=137
1	822565	rs111427246	C	T	.	.	ENZ=(SfiI|GGCCNNNN^NGGCC|GGCCAGCTTGGcC|822554|+);GNO;OTHERKG;RS=111427246;RSPOS=822565;SAO=0;SLO;SSR=0;VC=SNV;VP=0x050100000001000102000100;WGT=1;dbSNPBuildID=132
1	837753	rs187865648	G	A	.	.	CAF=[0.9995,0.0004591];COMMON=0;ENZ=(SfiI|GGCCNNNN^NGGCC|gGCCATTCTGGCC|837753|+);KGPROD;KGPhase1;RS=187865648;RSPOS=837753;SAO=0;SSR=0;VC=SNV;VP=0x050000000001000014000100;WGT=1;dbSNPBuildID=135
1	839873	rs192553893	C	T	.	.	CAF=[0.6768,0.3232];COMMON=1;ENZ=(NotI|GC^GGCCGC|GcGGCCGC|839872|+);KGPROD;KGPhase1;OTHERKG;RS=192553893;RSPOS=839873;SAO=0;SSR=0;VC=SNV;VP=0x050000000001100016000100;WGT=1;dbSNPBuildID=135
1	839911	rs76652930	C	T	.	.	ASP;ENZ=(NotI|GC^GGCCGC|GcGGCCGC|839910|+);GNO;OTHERKG;RS=76652930;RSPOS=839911;SAO=0;SSR=0;VC=SNV;VP=0x050000000005000102000100;WGT=1;dbSNPBuildID=131
1	839912	rs369394889	GGCC	G	.	.	ENZ=(NotI|GC^GGCCGC|GCggccGC|839910|+);OTHERKG;RS=369394889;RSPOS=839913;SAO=0;SSR=0;VC=DIV;VP=0x050000000001000002000200;WGT=1;dbSNPBuildID=138
1	839933	rs146045242	C	T	.	.	CAF=[0.8852,0.1148];COMMON=1;ENZ=(NotI|GC^GGCCGC|GcGGCCGC|839932|+);KGPROD;KGPhase1;OTHERKG;RS=146045242;RSPOS=839933;SAO=0;SSR=0;VC=SNV;VP=0x050000000001100016000100;WGT=1;dbSNPBuildID=134
1	840009	rs140080750	C	T	.	.	CAF=[0.9692,0.03076];COMMON=1;ENZ=(NotI|GC^GGCCGC|GcGGCCGC|840008|+);KGPROD;KGPhase1;OTHERKG;RS=140080750;RSPOS=840009;SAO=0;SSR=0;VC=SNV;VP=0x050000000001100016000100;WGT=1;dbSNPBuildID=134
```
 
