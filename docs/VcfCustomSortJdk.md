# VcfCustomSortJdk

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Sort a  VCF with dynamically-compiled java expressions


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfsortjdk  [options] Files

Usage: vcfsortjdk [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --body
      user's code is the whole body of the filter class, not just the 'apply' 
      method. 
      Default: false
    -e, --expression
      The java code expression.
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --limit
      limit output to N variant. -1 == unlimited
      Default: -1
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    --nocode
       Don't show the generated code
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    --saveCodeInDir
      Save the generated java code in the following directory
    -f, --script
      The java source code file.
    --tag
      Insert the sorting index with this INFO tag. Useful if a standard sort 
      on coordinate is applied after.
      Default: IDX
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit

```


## Keywords

 * vcf
 * filter
 * java
 * jdk
 * sort



## Creation Date

20260426

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcffilterjs/VcfCustomSortJdk.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcffilterjs/VcfCustomSortJdk.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcffilterjs/VcfCustomSortJdkTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcffilterjs/VcfCustomSortJdkTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfsortjdk** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)




## About the script


The user code is a piece of java code that will be inserted as the method apply, or as the body of a com.github.lindenb.jvarkit.tools.vcffilterjs.VcfCustomSortJdk.AbstractComparator class with implements `java.util.Comparator &lt;VariantContext&gt;`:

At the time of writing the documentation, the parent class AbstractComparator is defined as:

```java
public static class AbstractComparator
	extends com.github.lindenb.jvarkit.util.vcf.VcfTools
	implements Comparator<VariantContext>>
	{
	protected final VCFHeader header;
	protected AbstractComparator(final VCFHeader header) {
		super(header);
		this.header = header;
		}
	@Override
	public int compare(final VariantContext vc1 , final VariantContext vc2) {
		throw new IllegalStateException("apply(variant) for AbstractComparator is not implemented");
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
 9  public class VcfFilterJdkCustom123 extends com.github.lindenb.jvarkit.tools.vcffilterjs.VcfCustomSortJdk.AbstractComparator {
10    public VcfFilterJdkCustom123(final VCFHeader header) {
11    super(header);
12    }
13    @Override
	public int compare(final VariantContext vc1 , final VariantContext vc2) {
15     // user code starts here 
16     user's code is inserted here <===================
17     // user code ends here 
18     }
19  }
```

When the option `--body` is set : the user's code is the whole body (but the constructor) of the class




