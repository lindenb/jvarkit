# VcfGrantham

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

add grantham score from annotated VCF variant


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfgrantham  [options] Files

Usage: vcfgrantham [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * vcf
 * grantham



## Creation Date

20230503

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfgrantham/VcfGrantham.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfgrantham/VcfGrantham.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfgrantham** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java -jar dist/jeter.jar src/test/resources/rotavirus_rf.ann.vcf.gz | grep GRANT -m3
##INFO=<ID=GRANTHAM_SCORE,Number=1,Type=Integer,Description="Max grantham score for this variant">
RF01	970	.	A	C	48.67	.	AC=2;AN=10;ANN=C|missense_variant|MODERATE|Gene_18_3284|Gene_18_3284|transcript|AAA47319.1|protein_coding|1/1|c.952A>C|p.Lys318Gln|952/3267|952/3267|318/1088||;BQB=0.572843;DP=36;DP4=19,7,3,5;GRANTHAM_SCORE=53;HOB=0.32;ICB=0.425;MQ=60;MQ0F=0;MQB=1;MQSB=1;RPB=0.658863;SGB=10.3229;VDB=0.693968	GT:PL	0/0:0,9,47	0/0:0,18,73	0/0:0,18,73	0/0:0,33,116	1/1:95,24,0
RF02	251	.	A	T	21.29	.	AC=2;AN=10;ANN=T|stop_gained|HIGH|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32213.1|protein_coding|1/1|c.235A>T|p.Lys79*|235/2643|235/2643|79/880||,T|upstream_gene_variant|MODIFIER|Gene_1621_1636|Gene_1621_1636|transcript|CAA32214.1|protein_coding||c.-1371A>T|||||1371|WARNING_TRANSCRIPT_INCOMPLETE,T|upstream_gene_variant|MODIFIER|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32215.1|protein_coding||c.-1758A>T|||||1758|WARNING_TRANSCRIPT_NO_START_CODON;BQB=1;DP=24;DP4=18,0,6,0;GRANTHAM_SCORE=255;HOB=0.08;ICB=0.235294;MQ=60;MQ0F=0;MQB=1;RPB=0.566154;SGB=2.05141;VDB=0.0744703	GT:PL	0/0:0,15,57	0/1:31,0,5	0/1:31,0,5	0/0:0,9,42	0/0:0,24,69
```

