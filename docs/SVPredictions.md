# SVPredictions

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Basic Variant Effect prediction using gtf


## Usage

```
Usage: svpredictions [options] Files
  Options:
    --bnd
      Ignore the INFO/END attribute for SVTYPE=BND, so it is just considered 
      as a single point mutation.
      Default: false
    -F, --filter
      FILTER to set if variant failing prediction of option --where. Empty: no 
      FILTER, discard variant.
      Default: BAD_SV_PRED
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
  * -g, --gtf
      A GTF (General Transfer Format) file. See 
      https://www.ensembl.org/info/website/upload/gff.html .
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --max-genes
      don't print the genes names if their count exceed 'x'. '-1' = 
      ignore/unlimited 
      Default: 20
    -nti, --no-transcript-id
      don't print transcript id (reduce length of annotation)
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -r, --remove-attribute
      Do not print the annotations that don't contain the contraint for the 
      argument  --where
      Default: false
    --tag
      VCF info attribute
      Default: SVCSQ
    -u, --upstream
      Gene Upstream/Downstream length. A distance specified as a positive 
      integer.Commas are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb 
      Default: 5000
    --version
      print version and exit
    -w, --where
      where in gene should overlap the variant. Empty string: no limit/use all 
      possible annotations. Should be a comma/space/semicolon string with the 
      following items: 
      'intergenic|gene|transcript|intron|exon|utr|cds|downstream|upstream' 
      Default: <empty string>

```


## Keywords

 * vcf
 * annotation
 * prediction
 * sv


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew svpredictions
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190815

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfannot/SVPredictions.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfannot/SVPredictions.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfannot/SVPredictionsTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfannot/SVPredictionsTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **svpredictions** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
java -Xmx3g -Djava.io.tmpdir=. -jar dist/svpredictions.jar --max-genes 30  --gtf "human.gtf.gz" in.vcf >   out.vcf

more out.vcf
(...)
chr19	54672382	MantaBND:2392:1:3:1:0:0:0	G	[chr9:87618877[G	.	.	BND_PAIR_COUNT=7;CIPOS=-135,135;CLUSTER=CTX3378;IMPRECISE;MATEID=MantaBND:2392:1:3:1:0:0:1;PAIR_COUNT=7;SVCSQ=upstream_transcript_variant|ENSG00000167608|ENST00000416963|TMC4|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000494594|TMC4|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000468343|TMC4|protein_coding,exon|ENSG00000167608|ENST00000446291|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000453320|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000414665|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000437868|MBOAT7|protein_coding,intron|ENSG00000167608|ENST00000479750|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000494142|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000391754|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000465790|TMC4|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000495398|TMC4|protein_coding,exon|ENSG00000167608|ENST00000476013|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000474910|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000449249|MBOAT7|protein_coding,cds|ENSG00000167608|ENST00000376591|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000338624|MBOAT7|protein_coding,cds|ENSG00000167608|ENST00000301187|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000495968|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000497518|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000491216|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000245615|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000449860|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000495279|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000464098|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000431666|MBOAT7|protein_coding;SVTYPE=BND
chr21	10475514	MantaINS:141610:0:0:0:1:0	AG	AAAAAAAAAAAAAAA	.	.	CIGAR=1M14I1D;CLUSTER=CTX3514;DOWNSTREAM_PAIR_COUNT=0;END=10475515;PAIR_COUNT=0;SVCSQ=exon|ENSG00000270533|ENST00000604687|bP-21201H5.1|pseudogene;SVLEN=14;SVTYPE=INS;UPSTREAM_PAIR_COUNT=0
chr22	23478420	MantaDEL:144501:0:1:0:0:0	T	<DEL>	.	.	CIEND=-160,160;CIPOS=-174,175;CLUSTER=CTX3616;DOWNSTREAM_PAIR_COUNT=16;END=23479619;IMPRECISE;PAIR_COUNT=16;SVCSQ=utr&cds&intron&exon|ENSG00000100218|ENST00000406876|RTDR1|protein_coding,intron|ENSG00000100218|ENST00000216036|RTDR1|protein_coding,upstream_transcript_variant|ENSG00000272019|ENST00000606537|Metazoa_SRP|misc_RNA,transcript_ablation|ENSG00000221069|ENST00000408142|AC000029.1|miRNA,intron|ENSG00000100218|ENST00000439064|RTDR1|protein_coding,upstream_transcript_variant|ENSG00000100218|ENST00000421213|RTDR1|protein_coding,utr&intron&exon|ENSG00000100218|ENST00000452757|RTDR1|protein_coding;SVLEN=-1199;SVTYPE=DEL;UPSTREAM_PAIR_COUNT=16
(...)

```


