# VcfScanUpstreamOrf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Scan BAM for upstream-ORF. Inspired from https://github.com/ImperialCardioGenetics/uORFs 


## Usage

```
Usage: vcfscanupstreamorf [options] Files
  Options:
    --canonical
      reduce the number of transcripts. Keep one if some share the same UTR
      Default: false
    --dac
      disable scan for ATG creation
      Default: false
    --dad
      disable scan for ATG deletion
      Default: false
    --dsc
      disable scan for STOP creation
      Default: false
    --dsd
      disable scan for STOP deletion
      Default: false
    --exclude-cds
      remove a uORF it if enterely overlaps a coding region of the exon of an 
      alternative transcript.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --kal
      disable scan for Kozak change
      Default: false
    -k, -K, --kg, -kg
      UCSC knownGene File/URL. The knowGene format is a compact alternative to 
      GFF/GTF because one transcript is described using only one line.	Beware 
      chromosome names are formatted the same as your REFERENCE. A typical 
      KnownGene file is 
      http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz 
      .If you only have a gff file, you can try to generate a knownGene file 
      with [http://lindenb.github.io/jvarkit/Gff2KnownGene.html](http://lindenb.github.io/jvarkit/Gff2KnownGene.html)
      Default: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeBasicV19.txt.gz
    -o, --output
      Output file. Optional . Default: stdout. If there is no argument and 
      'output' is a directory or it ends with '.zip' an archive containing the 
      fasta+bed of the uORF is created and the program exits.
  * -r, -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --strong
      only accept events that are related to 'Strong' Kozak pattern.
      Default: false
    --uorf-only
      only print variants having something to say about an uorf
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * uorf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfscanupstreamorf
```

The java jar file will be installed in the `dist` directory.


## Creation Date

2019-02-18

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/upstreamorf/VcfScanUpstreamOrf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/upstreamorf/VcfScanUpstreamOrf.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/upstreamorf/VcfScanUpstreamOrfTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/upstreamorf/VcfScanUpstreamOrfTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfscanupstreamorf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## inspiration

part of this code was inspired from: https://github.com/ImperialCardioGenetics/uORFs/blob/master/5primeUTRannotator/five_prime_UTR_annotator.pm

## Examples

### Example 1

```
 wget -q -O - "https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr1.vcf.bgz" |\
 bcftools annotate -x "INFO,FILTER" |\
 java -jar /home/lindenb/src/jvarkit-git/dist/vcfscanupstreamorf.jar \
 	-R human_g1k_v37.fasta  --uorf-only  --canonical  

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	89333	rs1008713359	A	G	283.15	.	UORF_DEL_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:C|atg-pos:89334|cap-atg:2708|atg-cds:40|atg-frame:atg-out-of-cds-frame|kozak-seq:GAAATGC|kozak-strength:Moderate|stop-frame:not-in-frame-stop|stop-pos:89318|atg-stop:16|pep:.
1	89359	rs1327179626	C	T	3839.47	.	UORF_ADD_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:A|atg-pos:89359|cap-atg:2683|atg-cds:65|atg-frame:atg-out-of-cds-frame|kozak-seq:TGCATGT|kozak-strength:Weak|stop-frame:in-frame-stop|stop-pos:89356|atg-stop:3|pep:M
1	89391	rs1332733110	T	C	2045.80	.	UORF_DEL_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:G|atg-pos:89391|cap-atg:2651|atg-cds:97|atg-frame:atg-out-of-cds-frame|kozak-seq:TGAATGA|kozak-strength:Weak|stop-frame:in-frame-stop|stop-pos:89382|atg-stop:9|pep:VNK
1	89555	rs1200434471	A	C	283.62	.	UORF_ADD_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:G|atg-pos:89557|cap-atg:2485|atg-cds:263|atg-frame:atg-out-of-cds-frame|kozak-seq:GAAATGA|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89452|atg-stop:105|pep:MKSQNVSQKIIYNVCVRKRQYPSNFESLHQKENSK,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:G|atg-pos:89557|cap-atg:1312|atg-cds:7|atg-frame:atg-out-of-cds-frame|kozak-seq:GAAATGA|kozak-strength:Moderate|stop-frame:not-in-frame-stop|stop-pos:89556|atg-stop:1|pep:.
1	89560	rs1234719556	C	A	448.62	.	UORF_DEL_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:T|atg-pos:89562|cap-atg:2480|atg-cds:268|atg-frame:atg-out-of-cds-frame|kozak-seq:TGAATGA|kozak-strength:Weak|stop-frame:in-frame-stop|stop-pos:89541|atg-stop:21|pep:IKLKVKM,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:T|atg-pos:89562|cap-atg:1307|atg-cds:12|atg-frame:atg-in-cds-frame|kozak-seq:TGAATGA|kozak-strength:Weak|stop-frame:not-in-frame-stop|stop-pos:89555|atg-stop:7|pep:.
1	89624	rs1166058274	T	C	606.05	.	UORF_DEL_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:G|atg-pos:89624|cap-atg:2418|atg-cds:330|atg-frame:atg-in-cds-frame|kozak-seq:ACAATGA|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89609|atg-stop:15|pep:VKELF,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:G|atg-pos:89624|cap-atg:1245|atg-cds:74|atg-frame:atg-out-of-cds-frame|kozak-seq:ACAATGA|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89609|atg-stop:15|pep:VKELF
1	89718	rs865856422	A	G	1466.11	.	UORF_DEL_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:C|atg-pos:89719|cap-atg:2323|atg-cds:425|atg-frame:atg-out-of-cds-frame|kozak-seq:AAAATGA|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89710|atg-stop:9|pep:TKL,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:C|atg-pos:89719|cap-atg:1150|atg-cds:169|atg-frame:atg-out-of-cds-frame|kozak-seq:AAAATGA|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89710|atg-stop:9|pep:TKL
1	89831	rs1209426147	A	G	372.62	.	UORF_DEL_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:C|atg-pos:89832|cap-atg:2210|atg-cds:538|atg-frame:atg-out-of-cds-frame|kozak-seq:CTTATGT|kozak-strength:Weak|stop-frame:in-frame-stop|stop-pos:89811|atg-stop:21|pep:TFAIYHT,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:C|atg-pos:89832|cap-atg:1037|atg-cds:282|atg-frame:atg-in-cds-frame|kozak-seq:CTTATGT|kozak-strength:Weak|stop-frame:in-frame-stop|stop-pos:89811|atg-stop:21|pep:TFAIYHT
1	89945	rs1376722481	G	C	297.51	.	UORF_ADD_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:G|atg-pos:89947|cap-atg:2095|atg-cds:653|atg-frame:atg-out-of-cds-frame|kozak-seq:AATATGC|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89803|atg-stop:144|pep:MPLASVSHLAKPRLRSGKMEAISSWERRQRRWEYYVATYVCNLPYLAL,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:G|atg-pos:89947|cap-atg:922|atg-cds:397|atg-frame:atg-out-of-cds-frame|kozak-seq:AATATGC|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89803|atg-stop:144|pep:MPLASVSHLAKPRLRSGKMEAISSWERRQRRWEYYVATYVCNLPYLAL
1	90032	rs866094671	C	T	14378.50	.	UORF_ADD_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:A|atg-pos:90032|cap-atg:2010|atg-cds:738|atg-frame:atg-in-cds-frame|kozak-seq:TTCATGG|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89951|atg-stop:81|pep:MGQLVSRAARETKPQCTFYSLCAHQTC,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:A|atg-pos:90032|cap-atg:837|atg-cds:482|atg-frame:atg-out-of-cds-frame|kozak-seq:TTCATGG|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89951|atg-stop:81|pep:MGQLVSRAARETKPQCTFYSLCAHQTC
```

### extract to bed format

```
track name="uORF" description="uORF for http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeBasicV19.txt.gz"
#chrom	chromStart	chromEnd	name	score	strand	thickStart	thickEnd	itemRgb	blockCount	blockSizes	blockStarts
chr1	34553	36081	ENST00000417324.1.uorf	100	-	35140	35736	0,0,255	1	1528	0
chr1	89294	120932	ENST00000466430.1.uorf	500	-	91254	91491	0,255,0	1	31638	0
chr1	89550	91105	ENST00000495576.1.uorf	100	-	90431	90590	0,0,255	1	1555	0
chr1	139789	140339	ENST00000493797.1.uorf	500	-	139816	140223	0,255,0	1	550	0
chr1	141473	149707	ENST00000484859.1.uorf	1000	-	146708	146978	255,0,0	1	8234	0
chr1	142807	146831	ENST00000490997.1.uorf	1000	-	142988	146482	255,0,0	1	4024	0
chr1	157783	157887	ENST00000410691.1.uorf	0	-	157848	157887	0,0,0	1	104	0
chr1	236111	267253	ENST00000424587.2.uorf	100	-	236759	236918	0,0,255	1	31142	0
chr1	453632	460480	ENST00000450983.1.uorf	1000	-	453980	454166	255,0,0	1	6848	0
chr1	521368	523833	ENST00000417636.1.uorf	500	-	522285	523620	0,255,0	1	2465	0
chr1	529838	532878	ENST00000357876.5.uorf	500	-	530001	532684	0,255,0	1	3040	0
chr1	562756	564390	ENST00000452176.1.uorf	500	-	562878	562995	0,255,0	1	1634	0
chr1	646721	655580	ENST00000414688.1.uorf	500	-	647189	655553	0,255,0	1	8859	0
chr1	677192	685396	ENST00000416385.1.uorf	100	-	682910	683180	0,0,255	1	8204	0
chr1	693612	693716	ENST00000411249.1.uorf	0	-	693689	693716	0,0,0	1	104	0
chr1	694411	700305	ENST00000417659.1.uorf	100	-	700133	700208	0,0,255	1	5894	0
chr1	700236	714006	ENST00000428504.1.uorf	500	-	705034	709660	0,255,0	1	13770	0
chr1	736258	745541	ENST00000447500.1.uorf	500	-	741231	745515	0,255,0	1	9283	0
chr1	745488	753092	ENST00000435300.1.uorf	500	-	752900	753047	0,255,0	1	7604	0
chr1	761585	762902	ENST00000473798.1.uorf	500	-	762082	762571	0,255,0	1	1317	0
chr1	803450	812283	ENST00000446136.1.uorf	100	-	810390	812268	0,0,255	1	8833	0
chr1	852249	855072	ENST00000417705.1.uorf	100	-	852976	854794	0,0,255	1	2823	0
chr1	889805	894689	ENST00000487214.1.uorf	1000	-	889839	894620	255,0,0	1	4884	0
chr1	916546	917473	ENST00000341290.2.uorf	1000	-	916549	917473	255,0,0	1	927	0
chr1	931345	933431	ENST00000606034.1.uorf	100	-	931510	932137	0,0,255	1	2086	0
chr1	935353	935552	ENST00000428771.2.uorf	500	-	935487	935544	0,255,0	1	199	0
chr1	947376	948573	ENST00000458555.1.uorf	100	-	947459	947507	0,0,255	1	1197	0
chr1	997587	998668	ENST00000442292.2.uorf	500	-	997810	998119	0,255,0	1	1081	0
chr1	1019305	1051623	ENST00000482816.1.uorf	1000	-	1019401	1026923	255,0,0	1	32318	0
chr1	1026923	1027554	ENST00000379320.1.uorf	100	-	1027028	1027400	0,0,255	1	631	0
chr1	1026923	1041507	ENST00000379319.1.uorf	100	-	1041338	1041410	0,0,255	1	14584	0
(...)
```

