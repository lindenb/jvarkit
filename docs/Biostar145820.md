# Biostar145820

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

subsample/shuffle BAM to fixed number of alignments.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar biostar145820  [options] Files

Usage: biostar145820 [options] Files
  Options:
    --bamcompression
      Compression Level. 0: no compression. 9: max compression;
      Default: 5
    -f, --filter, --jexl
      A JEXL Expression that will be used to filter out some sam-records (see 
      https://software.broadinstitute.org/gatk/documentation/article.php?id=1255). 
      An expression should return a boolean value (true=exclude, false=keep 
      the read). An empty expression keeps everything. The variable 'record' 
      is the current observed read, an instance of SAMRecord (https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html).
      Default: 'Accept all' (Empty expression)
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --output
      Output file. Optional . Default: stdout
    --reference, -R
      For reading/writing CRAM files. Indexed fasta Reference file. This file 
      must be indexed with samtools faidx and with picard/gatk 
      CreateSequenceDictionary or samtools dict
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --seed
      Random seed. -1 == current time
      Default: -1
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit
    -n
       number of reads. negative: all reads, shuffle output.
      Default: -1

```


## Keywords

 * sam
 * bam
 * shuffle



## See also in Biostars

 * [https://www.biostars.org/p/145820](https://www.biostars.org/p/145820)



## Creation Date

20150615

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar145820.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar145820.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar145820Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar145820Test.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar145820** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

 * MED25 connects enhancer-promoter looping and MYC2-dependent activation of jasmonate signalling. Wang et al. Nature Plants 5, 616-625 (2019)  https://doi.org/10.1038/s41477-019-0441-9 


## Example

```bash
$ java -jar jvarkit.jar biostar145820 -n 10  -o out.bam  in.bam
```

## CIted In:

 * "MED25 connects enhancer-promoter looping and MYC2-dependent activation of jasmonate signalling. Wang et al. Nature Plants 5, 616-625 (2019)  https://doi.org/10.1038/s41477-019-0441-9 
 * Makamure, C.E., Justinen, S., Martínez, D.E. et al. Cool temperature inhibits binary fission and results in phenotypic and transcriptomic changes that suggest inducible aging in Diadumene lineata. BMC Res Notes 18, 293 (2025). https://doi.org/10.1186/s13104-025-07378-x 
```

