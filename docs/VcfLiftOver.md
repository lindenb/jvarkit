# VcfLiftOver

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Lift-over a VCF file


## Usage

```
Usage: vcfliftover [options] Files
  Options:
    --adaptivematch
      Use adapative liftover minmatch using the ratio between the min allele 
      size and the longest allele size
      Default: false
  * -f, --chain
      LiftOver file.
    --chainvalid
      Ignore LiftOver chain validation
      Default: false
    -check, --check
      Check variant allele sequence is the same on REF
      Default: false
    -x, --failed
      (file.vcf) write variants failing the liftOver here. Optional.
    -failtag, --failtag
      failed INFO tag
      Default: LIFTOVER_FAILED
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --indel, --indels
      do not LiftOver indels
      Default: false
    --info
      remove attribute from INFO on the fly
      Default: []
    -m, --minmatch
      lift over min-match.
      Default: 0.95
    -o, --output
      Output file. Optional . Default: stdout
  * -D, -R, -r, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    -T, --tag
      INFO tag
      Default: LIFTOVER
    --version
      print version and exit

```


## Keywords

 * vcf
 * liftover


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfliftover
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/liftover/VcfLiftOver.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/liftover/VcfLiftOver.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/liftover/VcfLiftOverTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/liftover/VcfLiftOverTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfliftover** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

