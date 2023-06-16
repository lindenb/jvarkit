# VcfSpringFilter

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Uses the java spring Framework to build complex vcf filters


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfspringfilter  [options] Files

Usage: vcfspringfilter [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
  * -c, --config
      Spring XML configuration file (  https://docs.spring.io/spring-framework/docs/4.2.x/spring-framework-reference/html/xsd-configuration.html 
      ) .
      Default: []
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -m, --main
      Main bean name
      Default: main
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * vcf
 * java
 * spring
 * framework



## Creation Date

20230526

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfspring/VcfSpringFilter.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfspring/VcfSpringFilter.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfspringfilter** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# Motivation

Use the spring framework ( https://docs.spring.io/spring-framework/docs/4.2.x/spring-framework-reference/html/xsd-configuration.html )  to load VCF filters implementing `com.github.lindenb.jvarkit.variant.VariantAnnotator`


# XML Configuration example

```xml
<?xml version="1.0" encoding="UTF-8"?>
<beans xmlns="http://www.springframework.org/schema/beans"
    xmlns:context="http://www.springframework.org/schema/context"
    xmlns:util="http://www.springframework.org/schema/util"
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:schemaLocation="
        http://www.springframework.org/schema/beans
        http://www.springframework.org/schema/beans/spring-beans.xsd
        http://www.springframework.org/schema/util
        https://www.springframework.org/schema/util/spring-util.xsd
        "
	>
    <util:list id="main">
    	<bean class=" com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate" >
    		<constructor-arg index="0" type="String" value="!vc.isSNP()"/>
    		 <property name="softFiltering" value="true"/>
    		 <property name="filterName" value="HELLO"/>
    	</bean>
    	<bean class=" com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate" >
    		<constructor-arg index="0" type="String" value="vc.homVarCount &gt; 1"/>
    		 <property name="softFiltering" value="true"/>
    		 <property name="filterName" value="HELLO2"/>
    	</bean>
    </util:list>
</beans>
```

## Example:

```
$ java -jar dist/jvarkit.jar vcfspringfilter -c config.xml  src/test/resources/spring-variant-annotators.01.xml
##fileformat=VCFv4.2
(...)
##FILTER=<ID=HELLO,Description="Filtered with JEXL expression(s) ( A Java EXpression Language (JEXL) expressions to filter the variants from a VCF. Empty string will accept all variants. Expression returning a TRUE will accept the variant. See https://gatk.broadinstitute.org/hc/en-us/articles/360035891011  ) : !vc.isSNP()">
##FILTER=<ID=HELLO2,Description="Filtered with JEXL expression(s) ( A Java EXpression Language (JEXL) expressions to filter the variants from a VCF. Empty string will accept all variants. Expression returning a TRUE will accept the variant. See https://gatk.broadinstitute.org/hc/en-us/articles/360035891011  ) : vc.homVarCount > 1">
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
(...)
RF03	2573	.	A	G	17.83	HELLO	AC=4;AN=10;BQB=0.974597;DP=9;DP4=0,5,0,4;HOB=0.48;ICB=0.117361;MQ=60;MQ0F=0;MQB=0.974597;RPB=0.487298;SGB=0.473945;VDB=0.0516381	GT:PL	0/0:0,3,17	1/1:31,6,0	1/1:31,6,0	0/0:0,3,17	0/0:0,9,42
RF04	887	.	A	G	5.31	HELLO;HELLO2	AC=1;AN=10;BQB=1;DP=48;DP4=16,28,3,1;HOB=0.02;ICB=0.0439024;MQ=60;MQ0F=0;MQB=1;MQSB=1;RPB=0.90467;SGB=3.91248;VDB=0.811811	GT:PL	0/1:40,0,28	0/0:0,24,98	0/0:0,24,98	0/0:0,33,120	0/0:0,42,134
RF08	992	.	G	C	70	HELLO	AC=4;AN=10;BQB=1;DP=33;DP4=0,21,0,12;HOB=0.48;ICB=0.117361;MQ=60;MQ0F=0;MQB=1;RPB=0.73431;SGB=7.42075;VDB=0.750182	GT:PL	0/0:0,21,66	1/1:62,18,0	1/1:62,18,0	0/0:0,15,57	0/0:0,27,72
(...)
```



