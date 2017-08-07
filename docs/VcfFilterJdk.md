# VcfFilterJdk

Filtering VCF with in-memory-compiled java expressions


## Usage

```
Usage: vcffilterjdk [options] Files
  Options:
    --body
      user's code is the whole body of the filter class, not just the 'apply' 
      method. 
      Default: false
    -e, --expression
       (js expression). Optional.
    -F, --filter
      If not empty, variants won't be discarded and this name will be used in 
      the FILTER column
      Default: <empty string>
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    --nocode
       Don't show the generated code
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --saveCodeInDir
      Save the generated java code in the following directory
    -f, --script
       (js file). Optional.
    --version
      print version and exit

```


## Keywords

 * vcf
 * filter
 * java
 * jdk



## See also in Biostars

 * [https://www.biostars.org/p/266201](https://www.biostars.org/p/266201)


## Compilation

### Requirements / Dependencies

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make vcffilterjdk
```

The *.jar libraries are not included in the main jar file, so you shouldn't move them (https://github.com/lindenb/jvarkit/issues/15#issuecomment-140099011 ).
The required libraries will be downloaded and installed in the `dist` directory.

### edit 'local.mk' (optional)

The a file **local.mk** can be created edited to override/add some definitions.

For example it can be used to set the HTTP proxy:

```
http.proxy.host=your.host.com
http.proxy.port=124567
```
## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcffilterjs/VcfFilterJdk.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcffilterjs/VcfFilterJdk.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcffilterjdk** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)




## About the script


The user code is a piece of java code that will be inserted as the method apply, or as the body of a com.github.lindenb.jvarkit.tools.vcffilterjs.VcfFilterJdk.AbstractFilter class with implements `java.util.function.Function<VariantContext,Object>`:

At the time of writing the documentation, the parent class AbstractFilter is defined as:

```java
public static class AbstractFilter
	extends com.github.lindenb.jvarkit.util.vcf.VcfTools
	implements Function<VariantContext,Object>
	{
	protected final Map<String,Object> userData = new HashMap<>();
	protected final VCFHeader header;
	protected AbstractFilter(final VCFHeader header) {
		super(header);
		this.header = header;
		}
	@Override
	public Object apply(final VariantContext variant) {
		throw new IllegalStateException("apply(variant) for AbstractFilter is not implemented");
		}
	}
```

where

* the base class VcfTools contains some utilities for parsing VEP/SNPEFF annotations and detecting Mendelian Violations. see [https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/util/vcf/VcfTools.java](https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/util/vcf/VcfTools.java).
* 'header' is the VCF header
* 'userData' is a placeHolder where the user is free to put things.

'userData' will be filled with the following properties:

* `<"first.variant",Boolean>` the current variant is the first variant in the VCF
* `<"last.variant",Boolean>` the current variant is the last variant in the VCF

If the user puts `<"STOP",Boolean.TRUE>` in `userData` the scanning of the VCF will be aborted without error.

The user code will be inserted in the following java code:


```
 1  import java.util.*;
 2  import java.util.stream.*;
 3  import java.util.function.*;
 4  import htsjdk.samtools.util.*;
 5  import htsjdk.variant.variantcontext.*;
 6  import htsjdk.variant.vcf.*;
 7  import javax.annotation.Generated;
 8  @Generated("VcfFilterJdk")
 9  public class VcfFilterJdkCustom123 extends com.github.lindenb.jvarkit.tools.vcffilterjs.VcfFilterJdk.AbstractFilter {
10    public VcfFilterJdkCustom123(final VCFHeader header) {
11    super(header);
12    }
13    @Override
14    public Object apply(final VariantContext variant) {
15     // user code starts here 
16     user's code is inserted here <===================
17     // user code ends here 
18     }
19  }
```

When the option `--body` is set : the user's code is the whole body (but the constructor) of the class



The program is then compiled in **memory**.

The method `apply` returns an object that can either:

* a boolean : true accept the variant, false reject the variant
* a [VariantContext](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html) to replace the current variant
* a [java.util.List](https://docs.oracle.com/javase/8/docs/api/java/util/List.html)<[VariantContext](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html) > to replace the current variant with a list of variants.

## See also

* VcfFilterJS . Slower, using javascript syntax (rhino engine)

## About Galaxy

At first, this tool is not safe for a public Galaxy server, because the javascript code can access the filesystem.
But you can use the JVM parameter

```
-J-Djava.security.manager
```

to prevent it to access the filesystem. See [http://stackoverflow.com/questions/40177810](http://stackoverflow.com/questions/40177810)


##  Examples


###  Example

>  I'm interested in finding all sites (regardless of genotype call heterozygous or homozygous) where at least one of the alternative alleles have an AD value (Allelic Depth) greater than 10,

see [https://bioinformatics.stackexchange.com/questions/974/](https://bioinformatics.stackexchange.com/questions/974/)

```
java -jar dist/vcffilterjdk.jar -e 'return variant.getGenotypes().stream().filter(G->G.hasAD() && java.util.Arrays.stream(G.getAD()).skip(1).filter(AD->AD>10).findAny().isPresent()).findAny().isPresent();' 
```

###  Example

filter homozygotes for sample NA12878


```
$ curl -sL "https://raw.github.com/jamescasbon/PyVCF/master/vcf/test/gatk.vcf" |\
	java -jar dist/vcffilterjdk.jar  -e 'return variant.getGenotype("NA12878").isHom();' | grep CHROM -A2

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	BLANK	NA12878	NA12891	NA12892	NA19238	NA19239	NA19240
chr22	42522755	.	C	G	36.98	.	AC=1;AF=0.071;AN=14;BaseQRankSum=-14.866;DP=1527;DS;Dels=0.01;FS=0.000;HRun=0;HaplotypeScore=253.4254;MQ=197.36;MQ0=2;MQRankSum=-10.810;QD=0.15;ReadPosRankSum=-17.244	GT:AD:DP:GQ:PL	0/0:26,1:27:51:0,51,570	0/0:208,40:248:99:0,236,4169	0/0:192,56:249:99:0,114,4292	0/1:179,66:245:75:75,0,3683	0/0:214,32:246:99:0,172,4235	0/0:200,49:249:61:0,61,4049	0/0:195,50:246:32:0,32,3757
chr22	42523003	rs116917064	A	G	7113.55	.	AC=8;AF=0.571;AN=14;BaseQRankSum=6.026;DB;DP=1433;DS;Dels=0.00;FS=0.000;HRun=1;HaplotypeScore=101.7894;MQ=182.04;MQ0=0;MQRankSum=-2.501;QD=4.96;ReadPosRankSum=8.294	GT:AD:DP:GQ:PL	0/1:10,2:12:1:1,0,257	1/1:9,173:183:99:2385,273,0	0/1:153,95:249:99:355,0,2355	0/1:140,110:250:99:1334,0,2242	0/1:164,85:249:99:1070,0,2279	0/1:160,90:250:99:1245,0,2300	0/1:156,81:238:99:724,0,2764
```

### Example

first and second genotype are not the same:

```
java -jar dist/vcffilterjdk.jar -e 'return !variant.getGenotype(0).sameGenotype(variant.getGenotype(1));' 
```



