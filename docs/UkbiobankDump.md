# UkbiobankDump

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Dump data https://afb.ukbiobank.ac.uk/ ukbiobank server. The server might not like too many requests. Use a your own risk. Doesn't work with jdk17 (?!)


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar ukbbdump  [options] Files

Usage: ukbbdump [options] Files
  Options:
    --api
      API base url
      Default: https://afb.ukbiobank.ac.uk/api
    --debug
      debug
      Default: false
  * -R, --dict
      A SAM Sequence dictionary source: it can be a *.dict file, a fasta file 
      indexed with 'picard CreateSequenceDictionary' or 'samtools dict', or 
      any hts file containing a dictionary (VCF, BAM, CRAM, intervals...)
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --region
      limit to that region. An interval as the following syntax : 
      "chrom:start-end". Some jvarkit programs also allow the following syntax 
      : "chrom:middle+extend"  or "chrom:start-end+extend" or 
      "chrom:start-end+extend-percent%".A program might use a Reference 
      sequence to fix the chromosome name (e.g: 1->chr1)
    --retry
      max-retry
      Default: 10
    --seconds
      wait 'x' seconds between each call to the API
      Default: 10
    --version
      print version and exit
    -o
      Output file. Optional . Default: stdout

```


## Keywords

 * ukbiobank



## Creation Date

20250425

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/ukbiobank/UkbiobankDump.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/ukbiobank/UkbiobankDump.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **ukbbdump** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





# WARNING

I cannot explain why this software doesn't run with jdk17 (I always get a http error 403) but runs with jdk23...



```
$ java -jar jvarkit.jar ukbbdump -R ref.dict --region "chr3:38000000-39000000" 
Chrom	Pos	rsID	Ref	Alt	nHomozygotes	HGVSp	maxImpact	maxConsequence	alleleCount	alleleNum	alleleFreq	geneSymbol
chr3	38000013	.	A	G	0	NM_001370264.1:c.672+574A>G,NM_001370265.1:c.573+574A>G,NM_001385038.1:c.1182+574A>G,NM_001385039.1:c.1182+574A>G,NM_015873.4:c.1182+574A>G	LOWEST	intron_variant	13	981062	0.0000132509464233657	VILL
chr3	38000014	.	G	A	0	NM_001370264.1:c.672+575G>A,NM_001370265.1:c.573+575G>A,NM_001385038.1:c.1182+575G>A,NM_001385039.1:c.1182+575G>A,NM_015873.4:c.1182+575G>A	LOWEST	intron_variant	3	981092	0.0000030578172077644096	VILL
chr3	38000015	.	G	A	0	NM_001370264.1:c.672+576G>A,NM_001370265.1:c.573+576G>A,NM_001385038.1:c.1182+576G>A,NM_001385039.1:c.1182+576G>A,NM_015873.4:c.1182+576G>A	LOWEST	intron_variant	3	981092	0.0000030578172077644096	VILL
chr3	38000018	rs1319217247	G	A	0	NM_001370264.1:c.672+579G>A,NM_001370265.1:c.573+579G>A,NM_001385038.1:c.1182+579G>A,NM_001385039.1:c.1182+579G>A,NM_015873.4:c.1182+579G>A	LOWEST	intron_variant	2	981088	0.0000020385531165400047	VILL
chr3	38000021	rs6806791	G	A	32	NM_001370264.1:c.672+582G>A,NM_001370265.1:c.573+582G>A,NM_001385038.1:c.1182+582G>A,NM_001385039.1:c.1182+582G>A,NM_015873.4:c.1182+582G>A	LOWEST	intron_variant	1342	981080	0.0013678802951848984	VILL
chr3	38000022	rs937502348	T	C	0	NM_001370264.1:c.672+583T>C,NM_001370265.1:c.573+583T>C,NM_001385038.1:c.1182+583T>C,NM_001385039.1:c.1182+583T>C,NM_015873.4:c.1182+583T>C	LOWEST	intron_variant	4	981070	0.0000040771810370309966	VILL
chr3	38000027	.	T	C	0	NM_001370264.1:c.672+588T>C,NM_001370265.1:c.573+588T>C,NM_001385038.1:c.1182+588T>C,NM_001385039.1:c.1182+588T>C,NM_015873.4:c.1182+588T>C	LOWEST	intron_variant	1	981072	0.0000010192931813363342	VILL
chr3	38000033	.	C	G	0	NM_001370264.1:c.672+594C>G,NM_001370265.1:c.573+594C>G,NM_001385038.1:c.1182+594C>G,NM_001385039.1:c.1182+594C>G,NM_015873.4:c.1182+594C>G	LOWEST	intron_variant	3	981056	0.000003057929414834627	VILL
```
##


