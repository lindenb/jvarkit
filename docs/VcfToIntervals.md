# VcfToIntervals

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

split a vcf to interval or bed for parallelization


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcf2intervals  [options] Files

Usage: vcf2intervals [options] Files
  Options:
    --bed, --bed-output
      force BED format as output. (Default is '.interval_list')
      Default: false
    -D, --distance
      min size of an interval (or use option -N). A distance specified as a 
      positive integer.Commas are removed. The following suffixes are 
      interpreted : b,bp,k,kb,m,mb,g,gb
      Default: -1
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --intervals, --bed-input
      Search for intervals for EACH record of the provided bed file. VCF path 
      must be provided and indexed.
    --min-distance
      extends the interval if the last variant is withing distance 'x' of the 
      next interval. Ignore if negative.A distance specified as a positive 
      integer.Commas are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb,g,gb 
      Default: -1
    -N, --variants, --n-variants
      number of variants per interval (or use option -D)
      Default: -1
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * vcf
 * bed
 * interval



## See also in Biostars

 * [https://www.biostars.org/p/9506628](https://www.biostars.org/p/9506628)
 * [https://www.biostars.org/p/9529137](https://www.biostars.org/p/9529137)



## Creation Date

20211112

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcf2intervals/VcfToIntervals.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcf2intervals/VcfToIntervals.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcf2intervals/VcfToIntervalsTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcf2intervals/VcfToIntervalsTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcf2intervals** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

input is a VCF file or a VCF stream.
input must be sorted on chrom/pos.

## Example

```
$ java -jar dist/vcf2intervals.jar -N 5 src/test/resources/rotavirus_rf.vcf.gz  --min-distance 1
@HD	VN:1.6	SO:coordinate
@SQ	SN:RF01	LN:3302
@SQ	SN:RF02	LN:2687
@SQ	SN:RF03	LN:2592
@SQ	SN:RF04	LN:2362
@SQ	SN:RF05	LN:1579
@SQ	SN:RF06	LN:1356
@SQ	SN:RF07	LN:1074
@SQ	SN:RF08	LN:1059
@SQ	SN:RF09	LN:1062
@SQ	SN:RF10	LN:751
@SQ	SN:RF11	LN:666
@RG	ID:S1	SM:S1
@RG	ID:S2	SM:S2
@RG	ID:S3	SM:S3
@RG	ID:S4	SM:S4
@RG	ID:S5	SM:S5
@CO	vcf2intervals. compilation:20211112182935 githash:9b2ab03 htsjdk:2.24.1 date:20211112183125. cmd:-N 5 src/test/resources/rotavirus_rf.vcf.gz --min-distance 1
RF01	970	970	1	1
RF02	251	1965	5	1715
RF03	1221	2150	5	930
RF03	2201	2573	3	373
RF04	887	1860	5	974
RF04	1900	1920	2	21
RF05	41	1297	5	1257
RF05	1339	1339	1	1
RF06	517	1132	5	616
RF07	98	952	4	855
RF08	926	992	2	67
RF09	294	414	3	121
RF10	46	175	3	130
RF11	74	79	1	6
```

```
$ java -jar dist/vcf2intervals.jar --distance 300 --min-distance 0 src/test/resources/rotavirus_rf.vcf.gz  
@HD	VN:1.6	SO:coordinate
@SQ	SN:RF01	LN:3302
@SQ	SN:RF02	LN:2687
@SQ	SN:RF03	LN:2592
@SQ	SN:RF04	LN:2362
@SQ	SN:RF05	LN:1579
@SQ	SN:RF06	LN:1356
@SQ	SN:RF07	LN:1074
@SQ	SN:RF08	LN:1059
@SQ	SN:RF09	LN:1062
@SQ	SN:RF10	LN:751
@SQ	SN:RF11	LN:666
@RG	ID:S1	SM:S1
@RG	ID:S2	SM:S2
@RG	ID:S3	SM:S3
@RG	ID:S4	SM:S4
@RG	ID:S5	SM:S5
@CO	vcf2intervals. compilation:20211112182935 githash:9b2ab03 htsjdk:2.24.1 date:20211112183310. cmd:--distance 300 --min-distance 0 src/test/resources/rotavirus_rf.vcf.gz
RF01	970	970	1	1
RF02	251	251	1	1
RF02	578	877	2	300
RF02	1726	1965	2	240
RF03	1221	1242	2	22
RF03	1688	1708	2	21
RF03	2150	2315	3	166
RF03	2573	2573	1	1
RF04	887	991	2	105
RF04	1241	1262	2	22
RF04	1857	1920	3	64
RF05	41	41	1	1
RF05	499	795	2	297
RF05	879	879	1	1
RF05	1297	1339	2	43
RF06	517	695	4	179
RF06	1129	1132	1	4
RF07	98	225	2	128
RF07	684	952	2	269
RF08	926	992	2	67
RF09	294	414	3	121
RF10	46	175	3	130
RF11	74	79	1	6
```


