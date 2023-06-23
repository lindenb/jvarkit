# SamJdk

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Filters a BAM using a java expression compiled in memory.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar samjdk  [options] Files

Usage: samjdk [options] Files
  Options:
    --bamcompression
      Compression Level. 0: no compression. 9: max compression;
      Default: 5
    --body
      user's code is the whole body of the filter class, not just the 'apply' 
      method. 
      Default: false
    -e, --expression
      java expression
    -X, --fail
      Save dicarded reads in that file
    -f, --file
      java file. Either option -e or -f is required.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -N, --limit
      limit to 'N' records (for debugging).
      Default: -1
    --nocode
       Don't show the generated code
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --pair
      [20171110] PAIR-MODE .The signature of java function is `public Object 
      apply(final List<SAMRecord> records)`. This function must return `true` 
      to accept the whole list, `false` to reject eveything, or another 
      `List<SAMRecord>`.Input MUST be sorted on query name using picard 
      SortSam (not `samtools sort` 
      https://github.com/samtools/hts-specs/issues/5 ).
      Default: false
    -R, --reference
      For reading CRAM. Indexed fasta Reference file. This file must be 
      indexed with samtools faidx and with picard CreateSequenceDictionary
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --saveCodeInDir
      Save the generated java code in the following directory
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * java
 * jdk
 * filter



## See also in Biostars

 * [https://www.biostars.org/p/270879](https://www.biostars.org/p/270879)
 * [https://www.biostars.org/p/274183](https://www.biostars.org/p/274183)
 * [https://www.biostars.org/p/278902](https://www.biostars.org/p/278902)
 * [https://www.biostars.org/p/279535](https://www.biostars.org/p/279535)
 * [https://www.biostars.org/p/283969](https://www.biostars.org/p/283969)
 * [https://www.biostars.org/p/286284](https://www.biostars.org/p/286284)
 * [https://www.biostars.org/p/286585](https://www.biostars.org/p/286585)
 * [https://www.biostars.org/p/286851](https://www.biostars.org/p/286851)
 * [https://www.biostars.org/p/286819](https://www.biostars.org/p/286819)
 * [https://www.biostars.org/p/287057](https://www.biostars.org/p/287057)
 * [https://www.biostars.org/p/299673](https://www.biostars.org/p/299673)
 * [https://www.biostars.org/p/301080](https://www.biostars.org/p/301080)
 * [https://www.biostars.org/p/305526](https://www.biostars.org/p/305526)
 * [https://www.biostars.org/p/306034](https://www.biostars.org/p/306034)
 * [https://www.biostars.org/p/309143](https://www.biostars.org/p/309143)
 * [https://www.biostars.org/p/327317](https://www.biostars.org/p/327317)
 * [https://www.biostars.org/p/335998](https://www.biostars.org/p/335998)
 * [https://www.biostars.org/p/336965](https://www.biostars.org/p/336965)
 * [https://www.biostars.org/p/340479](https://www.biostars.org/p/340479)
 * [https://www.biostars.org/p/342675](https://www.biostars.org/p/342675)
 * [https://www.biostars.org/p/345679](https://www.biostars.org/p/345679)
 * [https://www.biostars.org/p/362298](https://www.biostars.org/p/362298)
 * [https://www.biostars.org/p/368754](https://www.biostars.org/p/368754)
 * [https://www.biostars.org/p/378205](https://www.biostars.org/p/378205)
 * [https://www.biostars.org/p/408279](https://www.biostars.org/p/408279)
 * [https://www.biostars.org/p/417123](https://www.biostars.org/p/417123)
 * [https://www.biostars.org/p/427976](https://www.biostars.org/p/427976)
 * [https://www.biostars.org/p/424431](https://www.biostars.org/p/424431)
 * [https://www.biostars.org/p/450160](https://www.biostars.org/p/450160)
 * [https://www.biostars.org/p/9464312](https://www.biostars.org/p/9464312)
 * [https://www.biostars.org/p/9489815](https://www.biostars.org/p/9489815)
 * [https://www.biostars.org/p/9493510](https://www.biostars.org/p/9493510)
 * [https://www.biostars.org/p/9498170](https://www.biostars.org/p/9498170)
 * [https://www.biostars.org/p/9511167](https://www.biostars.org/p/9511167)
 * [https://www.biostars.org/p/9524098](https://www.biostars.org/p/9524098)
 * [https://www.biostars.org/p/9532167](https://www.biostars.org/p/9532167)
 * [https://www.biostars.org/p/9537698](https://www.biostars.org/p/9537698)
 * [https://www.biostars.org/p/9567318](https://www.biostars.org/p/9567318)



## Creation Date

20170807

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samjs/SamJdk.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samjs/SamJdk.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/samjs/SamJdkTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/samjs/SamJdkTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samjdk** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

 * "bioalcidae, samjs and vcffilterjs: object-oriented formatters and filters for bioinformatics files" . Bioinformatics, 2017. Pierre Lindenbaum & Richard Redon  [https://doi.org/10.1093/bioinformatics/btx734](https://doi.org/10.1093/bioinformatics/btx734).


