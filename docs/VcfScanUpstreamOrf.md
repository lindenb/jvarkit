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
      Output file. Optional . Default: stdout
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

## Example

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


