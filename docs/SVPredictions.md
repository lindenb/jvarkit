# SVPredictions

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Basic Variant Effect prediction using gtf


## Usage

```
Usage: svpredictions [options] Files
  Options:
  * -g, --gtf
      GTF File
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --max-genes
      don't print the genes names if their count exceed 'x'
      Default: 20
    -nti, --no-transcript-id
      don't print transcript id (reduce length of annotation)
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --tag
      VCF info attribute
      Default: SVCSQ
    -u, --upstream
      Upstream size. A distance specified as a positive integer.Comma are 
      removed. The following suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 5000
    --version
      print version and exit

```


## Keywords

 * vcf
 * annotation
 * prediction
 * sv


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew svpredictions
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190815

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfannot/SVPredictions.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfannot/SVPredictions.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **svpredictions** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



