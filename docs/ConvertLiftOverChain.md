# ConvertLiftOverChain

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert the contigs in a liftover chain to match another REFerence. (eg. to remove chr prefix, unknown chromosomes etc...)


## Usage

```
Usage: convertliftoverchain [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -R1, --ref1
      Source chain REFference mapping. Default : no conversion. Chromosome 
      mapping file. If the file looks like a NGS file (vcf, bam) the mapping 
      is extracted from a dictionary; Otherwise, it is interpreted as a 
      mapping file ( See https://github.com/dpryan79/ChromosomeMappings )
    -R2, --ref2
      Destination chain REFference mapping. Default : no conversion. Indexed 
      fasta Reference file. This file must be indexed with samtools faidx and 
      with picard CreateSequenceDictionary
    --version
      print version and exit

```


## Keywords

 * chain
 * liftover


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew convertliftoverchain
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190409

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/liftover/ConvertLiftOverChain.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/liftover/ConvertLiftOverChain.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/liftover/ConvertLiftOverChainTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/liftover/ConvertLiftOverChainTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **convertliftoverchain** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

