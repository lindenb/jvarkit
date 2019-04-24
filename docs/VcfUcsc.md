# VcfUcsc

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

annotate an VCF with mysql UCSC data


## Usage

```
Usage: vcfucsc [options] Files
  Options:
    -a, --accept
      JEXL expression used to convert a row to String. Empty=defaul/accept 
      all. See the manual.
      Default: <empty string>
    -D, --database
      mysql database name.
      Default: hg19
    -e, --expression
      JEXL expression used to convert a row to String. Empty=default. See the 
      manual. 
      Default: <empty string>
    -x, --extend
      Extend variant coordinates by 'x' bases.
      Default: 0
    -fi, --filterIn
      Set this FILTER if any item is found in the database
    -fo, --filterOut
      Set this FILTER if no item is found in the database
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -jdbc, --jdbc
      Java Database Connectivity (JDBC) URI
      Default: jdbc:mysql://genome-mysql.cse.ucsc.edu
    -o, --output
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


## JEXL expressions:

## Jexl expression

Filtering and conversion to string are performed using a JEXL expression. See https://commons.apache.org/proper/commons-jexl/reference/syntax.html

the following names are defined for the jexl context

 * **row** : a **ResultSet** https://docs.oracle.com/javase/8/docs/api/java/sql/ResultSet.html
 * **meta** : a **ResultSetMetaData** https://docs.oracle.com/javase/8/docs/api/java/sql/ResultSetMetaData.html
 * other fields are the names of the column in the table.



## History

20190424: switch to jexl expression
20180206: faster creating a prepared statement for each bin.size. fix chromContig

## Example


```
java -jar dist/vcfucsc.jar --table snp142 -e '${name}' input.vcf
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr3	124290753	.	G	C	579.77	.	UCSC_HG19_SNP142=rs145115089
chr3	124290943	.	A	G	491.77	.	UCSC_HG19_SNP142=rs7372055
chr3	124291069	.	G	A	266.77	.	UCSC_HG19_SNP142=rs7373767
chr3	124291171	.	C	CA	240.73	.	.
chr3	124291245	.	A	G	563.77	.	UCSC_HG19_SNP142=rs12695439
chr3	124291351	.	A	G	194.77	.	UCSC_HG19_SNP142=rs7613600
chr3	124291416	.	G	T	308.77	.	UCSC_HG19_SNP142=rs73189597
chr3	124291579	.	T	C	375.77	.	UCSC_HG19_SNP142=rs7649882
```
## Example

```
 java -jar dist/vcfucsc.jar --table vistaEnhancers  --tag VISTAENHANCERS -x 1000 -e '${chromStart}|${chromEnd}|${name}|${score}' input.vcf

```

