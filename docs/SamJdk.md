# SamJdk

Filters a BAM using a java expression compiled in memory.


## Usage

```
Usage: samjdk [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    --body
      user's code is the whole body of the filter class, not just the 'apply' 
      method. 
      Default: false
    -e, --expression
      javascript expression
    -X, --fail
      Save dicarded reads in that file
    -f, --file
      javascript file
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -N, --limit
      limit to 'N' records (for debugging).
      Default: -1
    --nocode
       Don't show the generated code
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --samoutputformat
      Sam output format.
      Default: TypeImpl{name='SAM', fileExtension='sam', indexExtension='null'}
    --saveCodeInDir
      Save the generated java code in the following directory
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * java
 * jdk
 * filter


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
$ make samjdk
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samjs/SamJdk.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samjs/SamJdk.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samjdk** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

Filters a BAM using a java expression compiled in memory.


## About the script


The user code is a piece of java code that will be inserted as the method test, or as the body of a com.github.lindenb.jvarkit.tools.samjs.SamJdk.AbstractFilter class with implements `java.util.function.Predicate<SAMRecord>`:

At the time of writing the documentation, the parent class AbstractFilter is defined as:

```java
public static class AbstractFilter
		implements Function<SAMRecord,Object>
		{
		// hashmap, the user is free to use It 
		protected final Map<String,Object> userData = new HashMap<>();
		// input SAM header
		protected final SAMFileHeader header;
		protected AbstractFilter(final SAMFileHeader header) {
			this.header = header;
			}
		@Override
		public Object apply(final SAMRecord record) {
			throw new IllegalStateException("apply(record) for AbstractFilter is not implemented");
			}
		}
```

where

* 'header' is the input SAM header
* 'userData' is a placeHolder where the user is free to put things.

The user code will be inserted in the following java code:


```
 1  import java.util.*;
 2  import java.util.stream.*;
 3  import java.util.function.*;
 4  import htsjdk.samtools.*;
 5  import htsjdk.samtools.util.*;
 6  import javax.annotation.Generated;
 7  @Generated(value="SamJdk",date="2017-08-07T14:48:39+0200")
 8  public class SamJdkCustom756098808 extends com.github.lindenb.jvarkit.tools.samjs.SamJdk.AbstractFilter {
 9    public SamJdkCustom756098808(final SAMFileHeader header) {
10    super(header);
11    }
12    @Override
13    public boolean test(final SAMRecord record) {
14     // user's code starts here 
15     return record.getContig()==null;
16     //user's code ends here 
17     }
18  }
```


When the option `--body` is set : the user's code is the whole body (but the constructor) of the class


The program is then compiled in **memory**.

The method `apply` returns an object that can either:

* a boolean : true accept the read, false reject the record
* a [SAMRecord](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html) to replace the current read
* a [java.util.List](https://docs.oracle.com/javase/8/docs/api/java/util/List.html)<[SAMRecord](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html) > to replace the current record with a list of records.


## Example

### Example 1


get a SAM where the  read OR the mate is unmapped

```bash
java -jar dist/samjdk.jar  \
	-e "return record.getReadUnmappedFlag() || record.getMateUnmappedFlag();" \
	ex1.bam

@HD	VN:1.4	SO:unsorted
@SQ	SN:seq1	LN:1575
@SQ	SN:seq2	LN:1584
B7_591:4:96:693:509	73	seq1	1	99	36M	*	0	0	CACTAGTGGCTCATTGTAAATGTGTGGTTTAACTCG	<<<<<<<<<<<<<<<;<<<<<<<<<5<<<<<;:<;7	H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:73
EAS54_65:7:152:368:113	73	seq1	3	99	35M	*	0	0	CTAGTGGCTCATTGTAAATGTGTGGTTTAACTCGT	<<<<<<<<<<0<<<<655<<7<<<:9<<3/:<6):H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:66
EAS51_64:8:5:734:57	137	seq1	5	99	35M	*	0	0	AGTGGCTCATTGTAAATGTGTGGTTTAACTCGTCC	<<<<<<<<<<<7;71<<;<;;<7;<<3;);3*8/5H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:66
B7_591:1:289:587:906	137	seq1	6	63	36M	*	0	0	GTGGCTCATTGTAATTTTTTGTTTTAACTCTTCTCT	(-&----,----)-)-),'--)---',+-,),''*,	H0:i:0	H1:i:0	MF:i:130	NM:i:5	UQ:i:38	Aq:i:63
EAS56_59:8:38:671:758	137	seq1	9	99	35M	*	0	0	GCTCATTGTAAATGTGTGGTTTAACTCGTCCATGG	<<<<<<<<<<<<<<<;<;7<<<<<<<<7<<;:<5%H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:72
EAS56_61:6:18:467:281	73	seq1	13	99	35M	*	0	0	ATTGTAAATGTGTGGTTTAACTCGTCCCTGGCCCA	<<<<<<<<;<<<8<<<<<;8:;6/686&;(16666H0:i:0	H1:i:1	MF:i:18	NM:i:1	UQ:i:5	Aq:i:39
EAS114_28:5:296:340:699	137	seq1	13	99	36M	*	0	0	ATTGTAAATGTGTGGTTTAACTCGTCCATGGCCCAG	<<<<<;<<<;<;<<<<<<<<<<<8<8<3<8;<;<0;	H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:73
B7_597:6:194:894:408	73	seq1	15	99	35M	*	0	0	TGTAAATGTGTGGTTTAACTCGTCCATTGCCCAGC	<<<<<<<<<7<<;<<<<;<<<7;;<<<*,;;572<H0:i:0	H1:i:1	MF:i:18	NM:i:1	UQ:i:9	Aq:i:43
EAS188_4:8:12:628:973	89	seq1	18	75	35M	*	0	0	AAATGTGTGGTTTAACTCGTCCATGGCCCAGCATT	==;=:;:;;:====;=;===:=======;==;===H0:i:1	H1:i:0	MF:i:64	NM:i:0	UQ:i:0	Aq:i:0
(...)
```

### Example 2

remove reads with indels:

```
java -jar dist/samjdk.jar -e 'if(record.getReadUnmappedFlag()) return false; Cigar cigar=record.getCigar();if(cigar==null) return false; for(int i=0;i< cigar.numCigarElements();++i) {if(cigar.getCigarElement(i).getOperator().isIndelOrSkippedRegion()) return false; } return true;' input.bam
```



