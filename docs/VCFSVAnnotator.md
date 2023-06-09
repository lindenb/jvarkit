# VCFSVAnnotator

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

SV Variant Effect prediction using gtf, gnomad, etc


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfsvannotator  [options] Files

Usage: vcfsvannotator [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -D, --define
      Dynamic parameters -Dkey=value .'extends': GTF Gene Upstream/Downstream 
      length. 'fraction': min common Fraction between two SVs/CNVs.
      Syntax: -Dkey=value
      Default: {extend=1000, fraction=0.9}
    --dgv
      DGV SV Variant file as Tabix indexed file from 
      http://dgv.tcag.ca/dgv/app/downloads 
    --ensemblreg
      GFF file from Ensembl, indexed with Tabix. e:g https://ftp.ensembl.org/pub/grch37/current/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff.gz
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    --gnomad
      Gnomad SV as BED file. Compressed with bgzip and indexed with tabix.
    --gtf
      GTF file
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    --regulomedb
      RegulomeDB bed sorted, bgzipped and indexed with tabix.
    --version
      print version and exit

```


## Keywords

 * vcf
 * annotation
 * prediction
 * sv
 * gnomad
 * gtf
 * regulome



## Creation Date

20190815

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfannot/VCFSVAnnotator.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfannot/VCFSVAnnotator.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfsvannotator** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
java -Xmx3g -Djava.io.tmpdir=. -jar dist/jvarkit.jar vcfsvannotator --gtf "human.gtf.gz" in.vcf >   out.vcf

more out.vcf
(...)
chr19	54672382	MantaBND:2392:1:3:1:0:0:0	G	[chr9:87618877[G	.	.	BND_PAIR_COUNT=7;CIPOS=-135,135;CLUSTER=CTX3378;IMPRECISE;MATEID=MantaBND:2392:1:3:1:0:0:1;PAIR_COUNT=7;SVCSQ=upstream_transcript_variant|ENSG00000167608|ENST00000416963|TMC4|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000494594|TMC4|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000468343|TMC4|protein_coding,exon|ENSG00000167608|ENST00000446291|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000453320|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000414665|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000437868|MBOAT7|protein_coding,intron|ENSG00000167608|ENST00000479750|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000494142|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000391754|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000465790|TMC4|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000495398|TMC4|protein_coding,exon|ENSG00000167608|ENST00000476013|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000474910|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000449249|MBOAT7|protein_coding,cds|ENSG00000167608|ENST00000376591|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000338624|MBOAT7|protein_coding,cds|ENSG00000167608|ENST00000301187|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000495968|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000497518|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000491216|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000245615|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000449860|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000495279|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000464098|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000431666|MBOAT7|protein_coding;SVTYPE=BND
chr21	10475514	MantaINS:141610:0:0:0:1:0	AG	AAAAAAAAAAAAAAA	.	.	CIGAR=1M14I1D;CLUSTER=CTX3514;DOWNSTREAM_PAIR_COUNT=0;END=10475515;PAIR_COUNT=0;SVCSQ=exon|ENSG00000270533|ENST00000604687|bP-21201H5.1|pseudogene;SVLEN=14;SVTYPE=INS;UPSTREAM_PAIR_COUNT=0
chr22	23478420	MantaDEL:144501:0:1:0:0:0	T	<DEL>	.	.	CIEND=-160,160;CIPOS=-174,175;CLUSTER=CTX3616;DOWNSTREAM_PAIR_COUNT=16;END=23479619;IMPRECISE;PAIR_COUNT=16;SVCSQ=utr&cds&intron&exon|ENSG00000100218|ENST00000406876|RTDR1|protein_coding,intron|ENSG00000100218|ENST00000216036|RTDR1|protein_coding,upstream_transcript_variant|ENSG00000272019|ENST00000606537|Metazoa_SRP|misc_RNA,transcript_ablation|ENSG00000221069|ENST00000408142|AC000029.1|miRNA,intron|ENSG00000100218|ENST00000439064|RTDR1|protein_coding,upstream_transcript_variant|ENSG00000100218|ENST00000421213|RTDR1|protein_coding,utr&intron&exon|ENSG00000100218|ENST00000452757|RTDR1|protein_coding;SVLEN=-1199;SVTYPE=DEL;UPSTREAM_PAIR_COUNT=16
(...)

```



