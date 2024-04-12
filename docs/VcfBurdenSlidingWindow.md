# VcfBurdenSlidingWindow

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

apply fisher test on VCF using a sliding window


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfburdenslidingwindow  [options] Files

Usage: vcfburdenslidingwindow [options] Files
  Options:
    --cases
      File or comma-separated list of control samples
    -C, --contig
      limit to this contig
    --controls
      File or comma-separated list of control samples
    -f, --filter
      A Java EXpression Language (JEXL) expressions to filter the variants 
      from a VCF. Empty string will accept all variants. Expression returning 
      a TRUE will accept the variant. See 
      https://gatk.broadinstitute.org/hc/en-us/articles/360035891011 
      Default: <empty string> (ACCEPT ALL)
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -save-vcf, --save-vcf
      Save Matching variants for the best p-value into that VCF.
    -t, --treshold
      fisher-test treshold. Discard results greater than this value.
      Default: 1.0
    --version
      print version and exit
    -s, --window-shift
      Window shift
      Default: 300
    -w, --window-size
      Window size
      Default: 1000

```


## Keywords

 * vcf
 * burden
 * case
 * control



## Creation Date

20190920

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenSlidingWindow.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenSlidingWindow.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenSlidingWindowTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenSlidingWindowTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfburdenslidingwindow** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# Motivation

apply fisher test on VCF using a sliding window

# Example

```
$ java -jar dist/jvarkit.jar vcfburdenslidingwindow --cases cases.txt --controls controls.txt -t 1 ./src/test/resources/test_vcf01.vcf  | head

#chrom	start0	end	name	length	p-value	affected_alt	affected_hom	unaffected_alt	unaffected_hom	variants.count
1	832199	833200	1:832200-833200	1001	1.0	0	3	1	2	1
1	832499	833500	1:832500-833500	1001	1.0	0	3	1	2	1
1	832799	833800	1:832800-833800	1001	1.0	0	3	1	2	1
1	839999	841000	1:840000-841000	1001	1.0	0	3	1	2	1
1	840299	841300	1:840300-841300	1001	1.0	0	3	1	2	1
1	840599	841600	1:840600-841600	1001	1.0	0	3	1	2	1
1	849299	850300	1:849300-850300	1001	1.0	0	3	1	2	1
1	849599	850600	1:849600-850600	1001	1.0	1	2	1	2	2
1	849899	850900	1:849900-850900	1001	1.0	1	2	1	2	2
```


