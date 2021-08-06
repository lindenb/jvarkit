# SwingBamCov

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Bam coverage viewer using Java Swing UI


## Usage

```
Usage: swingbamcov [options] Files
  Options:
    --gtf, --gff
      GFF3 file indexed with tabix to plot the genes.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -r, --regions, --interval
      default interval region on opening
      Default: <empty string>
    -q, --mapq
      min mapq
      Default: 1
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --small
      Display the reads when the region is small than 'x' bp. A distance 
      specified as a positive integer.Commas are removed. The following 
      suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: 200
    --version
      print version and exit

```


## Keywords

 * bam
 * alignment
 * graphics
 * visualization
 * swing


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew swingbamcov
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20210420

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfviewgui/SwingBamCov.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfviewgui/SwingBamCov.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **swingbamcov** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
java -jar dist/swingbamcov.jar -R ref.fa *.bam
```

## Screenshot

 * https://twitter.com/yokofakun/status/1392173415684100105

