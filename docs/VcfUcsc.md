# VcfUcsc

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

annotate an VCF with mysql UCSC data


## Usage

```
Usage: vcfucsc [options] Files
  Options:
    -a, --accept
      JEXL expression used to accept a result set. Must return a boolean. 
      Empty=defaul/accept all. See the manual. JEXL stands for Java EXpression 
      Language.  See 
      https://commons.apache.org/proper/commons-jexl/reference/syntax.html 
      Default: <empty string>
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -D, --database
      mysql database name.
      Default: hg19
    -e, --expression
      JEXL expression used to convert a row to String. Empty=default. See the 
      manual. JEXL stands for Java EXpression Language.  See 
      https://commons.apache.org/proper/commons-jexl/reference/syntax.html 
      Default: <empty string>
    -x, --extend
      Extend variant coordinates by 'x' bases.
      Default: 0
    -fi, --filterIn
      Set this FILTER if any item is found in the database
    -fo, --filterOut
      Set this FILTER if no item is found in the database
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -jdbc, --jdbc
      Java Database Connectivity (JDBC) URI
      Default: jdbc:mysql://genome-mysql.cse.ucsc.edu
    -L, --limit
      Limit number of items after filtration. Negative = no-limit
      Default: -1
    -o, --out
      Output file. Optional . Default: stdout
  * -T, -t, --table
      table name
    -tag, --tag
      INFO tag.
    --version
      print version and exit

```


## Keywords

 * ucsc
 * mysql
 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfucsc
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20160105

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfucsc/VcfUcsc.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfucsc/VcfUcsc.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfucsc/VcfUcscTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfucsc/VcfUcscTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfucsc** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

