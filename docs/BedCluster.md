# BedCluster

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Clusters a BED file into a set of BED files.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bedcluster  [options] Files

Usage: bedcluster [options] Files
  Options:
    -C, --contig, --chromosome
      group by chromosome.
      Default: false
    --consecutive
      When using option --size only use consecutive ordered regions. Default 
      is to find the best region anywhere.
      Default: false
    -F, --format
      An existing directory or a filename ending with the '.zip' or '.tar' or 
      '.tar.gz' suffix.
      Default: BED
      Possible Values: [BED, BED_GZ, INTERVAL_LIST, INTERVAL_LIST_GZ]
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -J, --jobs
      number of clusters. (or specify --size or --window-size/--window-shif)
      Default: -1
    -m, --manifest
      Manifest Bed file output containing chrom/start/end of each gene
    --md5-dir, --sub-dir
      prevent the creation of too many files in the same directory. Create 
      some intermediate directories based on filename's md5.
      Default: false
  * -o, --out
      An existing directory or a filename ending with the '.zip' or '.tar' or 
      '.tar.gz' suffix.
    -R, --reference
      For sorting and /or writing interval_list,A SAM Sequence dictionary 
      source: it can be a *.dict file, a fasta file indexed with 'picard 
      CreateSequenceDictionary' or 'samtools dict', or any hts file containing 
      a dictionary (VCF, BAM, CRAM, intervals...)
    -S, --size
      number of bases max per bin. (or specify --jobs or 
      --window-size/--window-shif). A distance specified as a positive 
      integer.Commas are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb,g,gb 
      Default: -1
    --version
      print version and exit
    --window-shift
      sliding  window shift (or use --size of --jobs).A distance specified as 
      a positive integer.Commas are removed. The following suffixes are 
      interpreted : b,bp,k,kb,m,mb,g,gb
      Default: -1
    --window-size
      sliding window size (or use --size of --jobs).A distance specified as a 
      positive integer.Commas are removed. The following suffixes are 
      interpreted : b,bp,k,kb,m,mb,g,gb
      Default: -1

```


## Keywords

 * bed
 * chromosome
 * contig



## See also in Biostars

 * [https://www.biostars.org/p/424828](https://www.biostars.org/p/424828)



## Creation Date

20200130

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bedcluster/BedCluster.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bedcluster/BedCluster.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bedcluster** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java -jar dist/jvarkit.jar bedcluster -j 10 -m jeter.mf -o jeter.zip --compress --contig test.bed

$ head jeter.mf

#chrom	start	end	filename	number_of_records	sum-length	avg-length	stddev-size
1	57460450	59012406	1_57460451_59012406.000000001.bed.gz	1	1551956	1551956	0
1	48998526	50489585	1_48998527_50489585.000000002.bed.gz	1	1491059	1491059	0
1	901876	248790491	1_901877_248790491.000000003.bed.gz	49	1393594	28440	41317
1	1567473	248154506	1_1567474_248154506.000000004.bed.gz	50	1393602	27872	39042
1	470970	229841608	1_470971_229841608.000000005.bed.gz	51	1393594	27325	37502
1	160445	248041507	1_160446_248041507.000000006.bed.gz	51	1393602	27325	37487
1	5647427	246685894	1_5647428_246685894.000000007.bed.gz	51	1393601	27325	37433
1	34553	245778447	1_34554_245778447.000000008.bed.gz	51	1393601	27325	37186
1	696290	247495148	1_696291_247495148.000000009.bed.gz	52	1393594	26799	35931

$ unzip -l jeter.zip | head
Archive:  jeter.zip
  Length      Date    Time    Name
---------  ---------- -----   ----
       76  2020-01-31 12:45   1_57460451_59012406.000000001.bed.gz
       74  2020-01-31 12:45   1_48998527_50489585.000000002.bed.gz
      487  2020-01-31 12:45   1_901877_248790491.000000003.bed.gz
      494  2020-01-31 12:45   1_1567474_248154506.000000004.bed.gz
      494  2020-01-31 12:45   1_470971_229841608.000000005.bed.gz
      511  2020-01-31 12:45   1_160446_248041507.000000006.bed.gz
      511  2020-01-31 12:45   1_5647428_246685894.000000007.bed.gz

 unzip -p jeter.zip 1_5647428_246685894.000000007.bed.gz | gunzip  -c  |head
1	5647427	5728355
1	8921060	8939308
1	11128527	11133154
1	13216355	13219581
1	16787442	16794976
1	20465818	20476879
1	21737952	21739786
1	26145130	26159432
1	37627163	37627235
1	38021842	38022108
```



