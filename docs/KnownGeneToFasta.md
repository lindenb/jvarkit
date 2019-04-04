# KnownGeneToFasta

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

convert knownGenes to fasta


## Usage

```
Usage: kg2fa [options] Files
  Options:
    --coding
      ignore non-coding transcripts.
      Default: false
    -D, --default
      Use default Known Gene source from UCSC.
      Default: false
    --dict
      Write optional dict file
    --empty
      Discard empty files.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --introns, --intron
      Remove introns
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
  * -r, -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --case, --style
      style: (0) do nothing (1): all UPPERCASE (2): all lowercase (3): exon 
      UPPERCASE + intron LOWERCASE . Otherwise do nothing
      Default: 0
    --utrs, --utr
      Remove UTRs
      Default: false
    --version
      print version and exit
    -L
      fasta line length.
      Default: 50

```


## Keywords

 * kg
 * knownGene
 * fasta


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew kg2fa
```

The java jar file will be installed in the `dist` directory.


## Creation Date

2019-02-13

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/kg2fa/KnownGeneToFasta.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/kg2fa/KnownGeneToFasta.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **kg2fa** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

