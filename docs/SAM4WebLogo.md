# SAM4WebLogo

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Sequence logo for different alleles or generated from SAM/BAM 


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar sam4weblogo  [options] Files

Usage: sam4weblogo [options] Files
  Options:
    -c, --clipped, --clip
      Use Clipped Bases
      Default: false
    --format, -F
      output format.
      Default: fasta
      Possible Values: [fasta, fastq, tabular]
    -fqp, --fqp
      [20180813] fastq padding quality character
      Default: -
    -fqu, --fqu
      [20180813] fastq unknown quality character
      Default: !
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -r, --region, --interval
      Region to observe. A source of intervals. The following suffixes are 
      recognized: vcf, vcf.gz bed, bed.gz, gtf, gff, gff.gz, gtf.gz.Otherwise 
      it could be an empty string (no interval) or a list of plain interval 
      separated by '[ \t\n;,]'
      Default: (unspecified)
    --naming
      How to print a Read. A format 'a la C-printf'. %% :% , %n:read name, %s: 
      read bases, %q: read quals, %f : read flags,%m: mapq, %c: contig, %b: 
      start, %B: unclipped start, %e: end, %E: unclipped end,%I: read group 
      id, %N: sample name,%S: SAM String.
      Default: %n (%f) %N
    --no-insert
      Do not show insertions
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -readFilter, --readFilter
      [20171201](moved to jexl)A JEXL Expression that will be used to filter 
      out some sam-records (see 
      https://software.broadinstitute.org/gatk/documentation/article.php?id=1255). 
      An expression should return a boolean value (true=exclude, false=keep 
      the read). An empty expression keeps everything. The variable 'record' 
      is the current observed read, an instance of SAMRecord (https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html).
      Default: 'Accept all' (Empty expression)
    -R, --reference
      For Reading CRAM. Indexed fasta Reference file. This file must be 
      indexed with samtools faidx and with picard/gatk 
      CreateSequenceDictionary or samtools dict
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * visualization
 * logo



## See also in Biostars

 * [https://www.biostars.org/p/73021](https://www.biostars.org/p/73021)
 * [https://www.biostars.org/p/368200](https://www.biostars.org/p/368200)
 * [https://www.biostars.org/p/9505463](https://www.biostars.org/p/9505463)
 * [https://www.biostars.org/p/9609826](https://www.biostars.org/p/9609826)
 * [https://www.biostars.org/p/103052](https://www.biostars.org/p/103052)



## Creation Date

20130524

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sam4weblogo/SAM4WebLogo.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sam4weblogo/SAM4WebLogo.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/sam4weblogo/SAM4WebLogoTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/sam4weblogo/SAM4WebLogoTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **sam4weblogo** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


