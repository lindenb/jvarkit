# SplitVcf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

split a vcf using a named list of intervals...


## Usage

```
Usage: splitvcf [options] Files
  Options:
    -g, --groupfile
      Chromosome group file. Intervals are 1-based. If undefined, splitvcf 
      will use the sequence dictionary to output one vcf per contig.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -m, --multi
      if set, allow one variant to be mapped on multiple chromosome group (the 
      record is duplicated)
      Default: false
  * -o, --out
      Output filename. Name must contain '__GROUPID__'
    -u, --unmapped
      unmapped interval name
      Default: OTHER
    --version
      print version and exit

```


## Keywords

 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew splitvcf
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SplitVcf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SplitVcf.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/SplitVcfTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/SplitVcfTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **splitvcf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ cat groups.txt

G1	10:112583204-112583210
G2	11
G3	12:1234-1235 13:20-30
```


```
$ java -jar dist/splitvcf.jar  -o tmp__GROUPID__.vcf.gz -g groups.txt in.vcf
$ ls tmp*
tmpG1.vcf.gz
tmpG2.vcf.gz
tmpG3.vcf.gz
tmpOTHER.vcf.gz
```


