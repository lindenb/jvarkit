# VCFBigWig

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Annotate a VCF with values from a bigwig file


## Usage

```
Usage: vcfbigwig [options] Files
  Options:
    -a, --aggregate
      How to aggregate overlapping values: 'avg' average; 'median': median, 
      'first': use first, 'all' : print all the data
      Default: avg
      Possible Values: [avg, median, first, all]
  * -B, --bigwig
      Path to the bigwig file. [20180122] If the path ends with '.xml' it is 
      interpretted as a XML file containing describing a set of BigWig 
      resources; See online doc.
    -C, --contained
      Specifies wig values must be contained by region. if false: return any 
      intersecting region values
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -T, --tag, -tag
      Name of the INFO tag. default: name of the bigwig
    --version
      print version and exit

```


## Keywords

 * vcf
 * wig
 * wiggle
 * bigwig


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfbigwig
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfbigwig/VCFBigWig.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfbigwig/VCFBigWig.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfbigwig/VCFBigWigTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfbigwig/VCFBigWigTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfbigwig** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## XML definition

multiple BigWig files can be specified using a XML file.

* Root is `<registry>`
* under `<registry>` is a set of `<bigwig>' elements.
* under `<bigwig>` contains the `<uri>'(required) , `<tag>` and `<description>`


## Example

```bash
 java -jar dist/vcfbigwig.jar \
 	-T GERP \
 	-B gerp.bw input.vcf.gz 
	
##INFO=<ID=GERP,Number=1,Type=Float,Description="Values from bigwig file: com.github.lindenb.jvarkit.tools.vcfbigwig.VCFBigWig BIGWIG=gerp.bw TAG=GERP IN=input.vcf.gz    VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO(...)
A	33926	.	G	A	182	.	GERP=-6.35(...)
A	45365	.	A	G	222	.	GERP=-3.55(...)
```


