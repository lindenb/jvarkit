# SamGrep

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

grep read-names in a bam file


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar samgrep  [options] Files

Usage: samgrep [options] Files
  Options:
    --bamcompression
      Compression Level. 0: no compression. 9: max compression;
      Default: 5
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -V, --invert
      invert
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    -f, --readfile
      file containing a list of read names
    -name, --readname
      add the read name
      Default: []
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --regions
      Limit analysis to this interval. A source of intervals. The following 
      suffixes are recognized: vcf, vcf.gz bed, bed.gz, gtf, gff, gff.gz, 
      gtf.gz.Otherwise it could be an empty string (no interval) or a list of 
      plain interval separated by '[ \t\n;,]'
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    -n, --stopafter
      when found, remove the read from the list of names when found more that 
      'n' time (increase speed)
      Default: -1
    --validation-stringency
      SAM Reader Validation Stringency
      Default: LENIENT
      Possible Values: [STRICT, LENIENT, SILENT]
    --version
      print version and exit

```


## Keywords

 * sam
 * bam



## Creation Date

20130506

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samgrep/SamGrep.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samgrep/SamGrep.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samgrep** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



### Examples


#### Example 1


```

java -jar  dist/samgrep.jar -R r001  -- samtools-0.1.18/examples/toy.sam 

@HD     VN:1.4  SO:unsorted
@SQ     SN:ref  LN:45
@SQ     SN:ref2 LN:40
@PG     ID:0    PN:com.github.lindenb.jvarkit.tools.samgrep.SamGrep     VN:dac03b80e9fd88a15648b22550e57d10c9bed725     CL:-R r001 samtools-0.1.18/examples/toy.sam
r001    163     ref     7       30      8M4I4M1D3M      =       37      39      TTAGATAAAGAGGATACTG     *       XX:B:S,12561,2,20,112
r001    83      ref     37      30      9M      =       7       -39     CAGCGCCAT       *

```


#### Example 4


```

java -jar  dist/samgrep.jar -R r001 -- -n 1 samtools-0.1.18/examples/toy.sam 

@HD     VN:1.4  SO:unsorted
@SQ     SN:ref  LN:45
@SQ     SN:ref2 LN:40
@PG     ID:0    PN:com.github.lindenb.jvarkit.tools.samgrep.SamGrep     VN:dac03b80e9fd88a15648b22550e57d10c9bed725     CL:-R r001 -n 1 samtools-0.1.18/examples/toy.sam
r001    163     ref     7       30      8M4I4M1D3M      =       37      39      TTAGATAAAGAGGATACTG     *       XX:B:S,12561,2,20,112

```







### Examples

#### Example 1

```

java -jar  dist/samgrep.jar -R r001  -- samtools-0.1.18/examples/toy.sam 

@HD	VN:1.4	SO:unsorted
@SQ	SN:ref	LN:45
@SQ	SN:ref2	LN:40
@PG	ID:0	PN:com.github.lindenb.jvarkit.tools.samgrep.SamGrep	VN:dac03b80e9fd88a15648b22550e57d10c9bed725	CL:-R r001 samtools-0.1.18/examples/toy.sam
r001	163	ref	7	30	8M4I4M1D3M	=	37	39	TTAGATAAAGAGGATACTG	*	XX:B:S,12561,2,20,112
r001	83	ref	37	30	9M	=	7	-39	CAGCGCCAT	*

```





#### Example 4



```

java -jar  dist/samgrep.jar -R r001 -- -n 1 samtools-0.1.18/examples/toy.sam 

@HD	VN:1.4	SO:unsorted
@SQ	SN:ref	LN:45
@SQ	SN:ref2	LN:40
@PG	ID:0	PN:com.github.lindenb.jvarkit.tools.samgrep.SamGrep	VN:dac03b80e9fd88a15648b22550e57d10c9bed725	CL:-R r001 -n 1 samtools-0.1.18/examples/toy.sam
r001	163	ref	7	30	8M4I4M1D3M	=	37	39	TTAGATAAAGAGGATACTG	*	XX:B:S,12561,2,20,112

```





