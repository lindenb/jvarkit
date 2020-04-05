# HtsFileServer

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Jetty Based http server serving Vcf and Bam files.


## Usage

```
Usage: htsfileserver [options] Files
  Options:
    --gtf
      Optional GTF file. Will be used to retrieve an interval by gene name
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -G, --no-genotype
      remove genotypes from vcf
      Default: false
    --port
      server port.
      Default: 8080
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --version
      print version and exit

```


## Keywords

 * vcf
 * bam
 * server



## See also in Biostars

 * [https://www.biostars.org/p/430718](https://www.biostars.org/p/430718)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew htsfileserver
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20200405

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/server/HtsFileServer.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/server/HtsFileServer.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **htsfileserver** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## input

input is a set of indexed Vcf/Bam file or a file with the suffix `.list` containing the path to the files.
 
## Example

```
java -jar dist/htsfileserver.jar -R src/test/resources/rotavirus_rf.fa src/test/resources/S*.bam  src/test/resources/rotavirus_rf.*.vcf.gz
```

# Screenshot

![https://i.imgur.com/ObRsVxE.png](https://i.imgur.com/ObRsVxE.png)

