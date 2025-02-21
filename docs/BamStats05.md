# BamStats05

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Coverage statistics for a BED file, group by gene


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bamstats05  [options] Files

Usage: bamstats05 [options] Files
  Options:
  * -B, --bed
      bed file (columns: chrom(tab)start(tab)end(tab)GENE)
    -f, --filter, --jexl
      A JEXL Expression that will be used to filter out some sam-records (see 
      https://software.broadinstitute.org/gatk/documentation/article.php?id=1255). 
      An expression should return a boolean value (true=exclude, false=keep 
      the read). An empty expression keeps everything. The variable 'record' 
      is the current observed read, an instance of SAMRecord (https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html).
      Default: record.getMappingQuality()<1 || record.getDuplicateReadFlag() || record.getReadFailsVendorQualityCheckFlag() || record.isSecondaryOrSupplementary()
    --groupby
      Group Reads by. Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --mapq
      Min mapping quality
      Default: 1
    -merge, --merge
      [20181122] Merge overlapping intervals for the same gene.
      Default: false
    -m, --mincoverage
      Coverage treshold. Any depth under this value will be considered as 
      'not-covered'.  Default: 0
      Default: []
    -o, --output
      Output file. Optional . Default: stdout
    -R, --reference
      For reading/writing CRAM files. Indexed fasta Reference file. This file 
      must be indexed with samtools faidx and with picard/gatk 
      CreateSequenceDictionary or samtools dict
    --version
      print version and exit

```


## Keywords

 * bam
 * coverage
 * statistics
 * bed



## See also in Biostars

 * [https://www.biostars.org/p/324639](https://www.biostars.org/p/324639)
 * [https://www.biostars.org/p/194393](https://www.biostars.org/p/194393)
 * [https://www.biostars.org/p/35083](https://www.biostars.org/p/35083)



## Creation Date

20151012

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bamstats04/BamStats05.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bamstats04/BamStats05.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bamstats04/BamStats05Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bamstats04/BamStats05Test.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bamstats05** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

input is one or more indexed bam file.

One file with  the suffix '.list' is interpreted as a text file with one path per line.

If there is no argument, stdin is interpreted as a list of path to the bam like in `find . -name "*.bam"`


## Cited In:

  * "Custom hereditary breast cancer gene panel selectively amplifies target genes for reliable variant calling" . BioRxiv https://doi.org/10.1101/322180
  * "Future MicrobiologyVol. 15, No. 15  Comparative genomics of Sporothrix species and identification of putative pathogenic-gene determinants.  Published Online:12 Nov 2020https://doi.org/10.2217/fmb-2019-0302"
  * "Key innovation triggers widespread radiation of the genus Medicago". Zhipeng Liu & al. https://doi.org/10.21203/rs.3.rs-3181566/v1
  * Over-expression and increased copy numbers of a cytochrome P450 and two UDP-glucuronosyltransferase genes in macrocyclic lactone resistant Psoroptes ovis of cattle Jack Hearn, Wouter van Mol, Roel Meyermans, Kathryn Bartley, Tyler Alioto, Jessica Gomez-Garrido, Fernando Cruz, Francisco Camara Ferreira, Marta Gut, Ivo Gut, Nadine Buys, Steven Janssens, Karyn Adams, Sara Roose,  Thomas Van Leeuwen, Wannes Dermauw, John S. Gilleard, Russell Avramenko, Peter Geldhof, Edwin Claerebout, Stewart T.G. Burgess bioRxiv 2025.02.06.636820; doi: https://doi.org/10.1101/2025.02.06.636820 


## Example

```
$ head genes.bed
1	179655424	179655582	ZORG
1	179656788	179656934	ZORG

$ java -jar  dist/bamstats05.jar -B genes.bed --mincoverage 10 in.bam > out.txt

$ head out.txt
#chrom	start	end	gene	sample	length	mincov	maxcov	avg	nocoverage.bp	percentcovered
1	179655424	179656934	ZORG	SAMPLE1	304	27	405	216.80921052631578	0	100
```



