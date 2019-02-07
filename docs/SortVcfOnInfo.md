# SortVcfOnInfo

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Sort a VCF a field in the INFO column


## Usage

```
Usage: sortvcfoninfo [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --output
      Output file. Optional . Default: stdout
  * -T, --tag, -t
      INFO tag
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit

```


## Keywords

 * vcf
 * sort
 * annotation


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew sortvcfoninfo
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sortvcfonref/SortVcfOnInfo.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sortvcfonref/SortVcfOnInfo.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **sortvcfoninfo** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


##Example

```bash
]$ curl  "https://raw.github.com/arq5x/gemini/master/test/test4.vep.snpeff.vcf" |\
   java -jar dist/sortvcfoninfo.jar -F BaseQRankSum | grep -vE "^#" 

chr10	1142208	.	T	C	3404.30	.	AC=8;AF=1.00;AN=8;
chr10	135336656	.	G	A	38.34	.	AC=4;AF=1.00;AN=4;
chr10	52004315	.	T	C	40.11	.	AC=4;AF=1.00;AN=4;
chr10	52497529	.	G	C	33.61	.	AC=4;AF=1.00;AN=4;
chr10	126678092	.	G	A	89.08	.	AC=1;AF=0.13;AN=8;BaseQRankSum=-3.120;
chr16	72057435	.	C	T	572.98	.	AC=1;AF=0.13;AN=8;BaseQRankSum=-2.270;
chr10	48003992	.	C	T	1047.87	.	AC=4;AF=0.50;AN=8;BaseQRankSum=-0.053;
chr10	135210791	.	T	C	65.41	.	AC=4;AF=0.50;AN=8;BaseQRankSum=2.054;
chr10	135369532	.	T	C	122.62	.	AC=2;AF=0.25;AN=8;BaseQRankSum=2.118;
```

