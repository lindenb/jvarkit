# FindCorruptedFiles

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Reads filename from stdin and prints corrupted NGS files (VCF/BAM/FASTQ/BED/TBI/BAI)


## Usage

```
Usage: findcorruptedfiles [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -E, --noempty
      empty file is an error.
      Default: false
    -V, --stringency
      BAM ValidationStringency
      Default: LENIENT
      Possible Values: [STRICT, LENIENT, SILENT]
    --version
      print version and exit
    -N
      number of features (samrecord, variant) to read. -1= read everything.
      Default: 100

```


## Keywords

 * vcf
 * bam
 * fastq
 * bed


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew findcorruptedfiles
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FindCorruptedFiles.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FindCorruptedFiles.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **findcorruptedfiles** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ find  DIR1 DIR2 -type f |\
java -jar dist/findcorruptedfiles.jar \
	-V SILENT 2> /dev/null > redo.txt
```

