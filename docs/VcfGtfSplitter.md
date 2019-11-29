# VcfGtfSplitter

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Split VCF+VEP by gene/transcript using a GTF file.


## Usage

```
Usage: vcfgtfsplitter [options] Files
  Options:
    --bcf
      Use bcf format
      Default: false
    -C, --contig, --chromosome
      Limit to those contigs.
      Default: []
    --coding
      Only use  gene_biotype="protein_coding".
      Default: false
    --upstream, --downstream
      length for upstream and downstream features. A distance specified as a 
      positive integer.Commas are removed. The following suffixes are 
      interpreted : b,bp,k,kb,m,mb
      Default: 1000
    --features
      Features to keep. Comma separated values. A set of 'cds,exon,intron,transcript,utr,utr5,utr3,stop,start,upstream,downstream,splice'
      Default: cds,exon,intron,transcript,cds_utr,cds_utr5,cds_utr3,utr5,utr3,stop,start
    --force
      Force writing a gene/transcript even if there is no variant.
      Default: false
  * -g, -G, --gtf
      A GTF (General Transfer Format) file. See 
      https://www.ensembl.org/info/website/upload/gff.html . Please note that 
      CDS are only detected if a start and stop codons are defined.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --ignore-filtered
      Ignore FILTERED variant
      Default: false
    --index
      index files
      Default: false
    -m, --manifest
      Manifest Bed file output containing chrom/start/end of each gene
  * -o, --output
      An existing directory or a filename ending with the '.zip' suffix.
    --splice
      distance to splice site for 'splice' feature. A distance specified as a 
      positive integer.Commas are removed. The following suffixes are 
      interpreted : b,bp,k,kb,m,mb
      Default: 5
    -T, --transcript
      split by transcript. (default is to split per gene)
      Default: false
    --version
      print version and exit
    --xannotate
      Remove annotations. Variant Attribute cleaner. The syntax is the same as 
      'bcftools annotate'. e.g: 'INFO/AC,INFO/ANN'  Empty string does nothing.

```


## Keywords

 * genes
 * vcf
 * split
 * gtf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfgtfsplitter
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20191118

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfGtfSplitter.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfGtfSplitter.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfGtfSplitterTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfGtfSplitterTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfgtfsplitter** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


### Example

```
$ java -jar dist/vcfgtfsplitter.jar  -m jeter.manifest --gtf  input.gtf.gz -o jeter.zip src/test/resources/test_vcf01.vcf 

$ unzip -l jeter.zip 
Archive:  jeter.zip
  Length      Date    Time    Name
---------  ---------- -----   ----
     4741  2019-11-18 10:52   d2/a0884d9bf86378a2cd0bfacef19723/ENSG00000188157.vcf.gz
     4332  2019-11-18 10:52   4d/11ff4f2413d4a369ee9a51192acad3/ENSG00000131591.vcf.gz
---------                     -------
     9073                     2 files

$ column -t jeter.manifest 
#chrom  start    end      Gene-Id          Gene-Name  Gene-Biotype    path                                                      Count_Variants
1       955502   991496   ENSG00000188157  AGRN       protein_coding  d2/a0884d9bf86378a2cd0bfacef19723/ENSG00000188157.vcf.gz  20
1       1017197  1051741  ENSG00000131591  C1orf159   protein_coding  4d/11ff4f2413d4a369ee9a51192acad3/ENSG00000131591.vcf.gz  6


$ java -jar dist/vcfgtfsplitter.jar -T --index  --gtf  jeter.gtf  -m jeter.manifest -o jeter.zip src/test/resources/test_vcf01.vcf
$ unzip -l jeter.zip 
Archive:  jeter.zip
  Length      Date    Time    Name
---------  ---------- -----   ----
     4749  2019-11-18 10:56   93/00d3aa560bdbe3b3ff964e24303107/ENST00000379370.vcf.gz
     4749  2019-11-18 10:56   93/00d3aa560bdbe3b3ff964e24303107/ENST00000379370.vcf.gz.tbi
     3108  2019-11-18 10:56   6d/001c37d11b192bcd77ac089ec258f2/ENST00000477585.vcf.gz
(...)
     2969  2019-11-18 10:56   3a/8e697b04e03716942b362afbe7ee51/ENST00000472741.vcf.gz
     2969  2019-11-18 10:56   3a/8e697b04e03716942b362afbe7ee51/ENST00000472741.vcf.gz.tbi
     2969  2019-11-18 10:56   2c/4ced556d7722b978453c38d5525ee0/ENST00000480643.vcf.gz
     2969  2019-11-18 10:56   2c/4ced556d7722b978453c38d5525ee0/ENST00000480643.vcf.gz.tbi
---------                     -------
   175030                     46 files

$ column -t jeter.manifest 
#chrom  start    end      Gene-Id          Gene-Name  Gene-Biotype    Transcript-Id    path                                                      Count_Variants
1       955502   991496   ENSG00000188157  AGRN       protein_coding  ENST00000379370  93/00d3aa560bdbe3b3ff964e24303107/ENST00000379370.vcf.gz  20
1       969485   976105   ENSG00000188157  AGRN       protein_coding  ENST00000477585  6d/001c37d11b192bcd77ac089ec258f2/ENST00000477585.vcf.gz  3
1       970246   976777   ENSG00000188157  AGRN       protein_coding  ENST00000469403  a9/fe6f84e1a425797d5b58ceae3a8313/ENST00000469403.vcf.gz  2
1       983908   984774   ENSG00000188157  AGRN       protein_coding  ENST00000492947  55/b5a514c94417efe2632aba379d8247/ENST00000492947.vcf.gz  1
1       1017197  1051461  ENSG00000131591  C1orf159   protein_coding  ENST00000379339  7b/26f9faff411f6e056fefc89cd15a24/ENST00000379339.vcf.gz  6
1       1017197  1051736  ENSG00000131591  C1orf159   protein_coding  ENST00000448924  56/fad8194b90771b0c3e931bc4772be1/ENST00000448924.vcf.gz  6
1       1017197  1051736  ENSG00000131591  C1orf159   protein_coding  ENST00000294576  df/3aa6bebac650a6c9f59409521ae17d/ENST00000294576.vcf.gz  6
(...)
```

# screenshot

* https://twitter.com/yokofakun/status/1197149666237911040

![https://twitter.com/yokofakun/status/1197149666237911040](https://pbs.twimg.com/media/EJ0hReMX0AcaBoq?format=png&name=small)

* https://twitter.com/yokofakun/status/1199621057533140992

![https://twitter.com/yokofakun/status/1199621057533140992](https://twitter.com/i/status/1199621057533140992)

