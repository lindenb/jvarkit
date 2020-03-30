# FastqSplit

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Split Fastq files into multiple files.


## Usage

```
Usage: fastqsplit [options] Files
  Options:
    -async, --async
      use async I/O
      Default: false
    -n, -N, --count
      number of reads per file. will be generated. Or use option 's'.
      Default: -1
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -ii, --input-interleaved
      input is paired interleaved
      Default: false
    -m, --manifest
      Optional manifest file
    -md5, --md5
      write md5 file
      Default: false
  * -o, --output
      Output file name. It MUST  end with a fastq suffix and MUST contain the 
      word __TOKEN__
    -oi, --output-interleaved
      output is interleaved
      Default: false
    -s, -S, --splits
      number of splits. At most 'x' (pair-of)files will be generated. Or use 
      option 'n'
      Default: -1
    -validate, --validate
      validate read names
      Default: false
    --version
      print version and exit

```


## Keywords

 * fastq


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew fastqsplit
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20200327

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/fastq/FastqSplit.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/fastq/FastqSplit.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/fastq/FastqSplitTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/fastq/FastqSplitTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **fastqsplit** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash


```

