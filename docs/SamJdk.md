# SamJdk

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

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
      java expression
    -X, --fail
      Save dicarded reads in that file
    -f, --file
      java file. Either option -e or -f is required.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -N, --limit
      limit to 'N' records (for debugging).
      Default: -1
    --nocode
       Don't show the generated code
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --pair
      [20171110] PAIR-MODE .The signature of java function is `public Object 
      apply(final List<SAMRecord> records)`. This function must return `true` 
      to accept the whole list, `false` to reject eveything, or another 
      `List<SAMRecord>`.Input MUST be sorted on query name using picard 
      SortSam (not `samtools sort` 
      https://github.com/samtools/hts-specs/issues/5 ).
      Default: false
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
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



## See also in Biostars

 * [https://www.biostars.org/p/270879](https://www.biostars.org/p/270879)
 * [https://www.biostars.org/p/274183](https://www.biostars.org/p/274183)
 * [https://www.biostars.org/p/278902](https://www.biostars.org/p/278902)
 * [https://www.biostars.org/p/279535](https://www.biostars.org/p/279535)
 * [https://www.biostars.org/p/283969](https://www.biostars.org/p/283969)
 * [https://www.biostars.org/p/286284](https://www.biostars.org/p/286284)
 * [https://www.biostars.org/p/286585](https://www.biostars.org/p/286585)
 * [https://www.biostars.org/p/286851](https://www.biostars.org/p/286851)
 * [https://www.biostars.org/p/286819](https://www.biostars.org/p/286819)
 * [https://www.biostars.org/p/287057](https://www.biostars.org/p/287057)
 * [https://www.biostars.org/p/299673](https://www.biostars.org/p/299673)
 * [https://www.biostars.org/p/301080](https://www.biostars.org/p/301080)
 * [https://www.biostars.org/p/305526](https://www.biostars.org/p/305526)
 * [https://www.biostars.org/p/306034](https://www.biostars.org/p/306034)
 * [https://www.biostars.org/p/309143](https://www.biostars.org/p/309143)
 * [https://www.biostars.org/p/327317](https://www.biostars.org/p/327317)
 * [https://www.biostars.org/p/335998](https://www.biostars.org/p/335998)
 * [https://www.biostars.org/p/336965](https://www.biostars.org/p/336965)
 * [https://www.biostars.org/p/340479](https://www.biostars.org/p/340479)
 * [https://www.biostars.org/p/342675](https://www.biostars.org/p/342675)
 * [https://www.biostars.org/p/345679](https://www.biostars.org/p/345679)
 * [https://www.biostars.org/p/362298](https://www.biostars.org/p/362298)
 * [https://www.biostars.org/p/368754](https://www.biostars.org/p/368754)
 * [https://www.biostars.org/p/378205](https://www.biostars.org/p/378205)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew samjdk
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samjs/SamJdk.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samjs/SamJdk.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/samjs/SamJdkTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/samjs/SamJdkTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samjdk** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

 * "bioalcidae, samjs and vcffilterjs: object-oriented formatters and filters for bioinformatics files" . Bioinformatics, 2017. Pierre Lindenbaum & Richard Redon  [https://doi.org/10.1093/bioinformatics/btx734](https://doi.org/10.1093/bioinformatics/btx734).


## Motivation

Filters a BAM using a java expression compiled at runtime.


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
 6  import javax.annotation.processing.Generated;
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

### Example 3

select unmappeds read or clipped reads

```
java -jar dist/samjdk.jar -o out.bam -e 'return record.getReadUnmappedFlag() || record.getCigar().getCigarElements().stream().anyMatch(C->C.getOperator().isClipping());'  in.bam
```

### Example 4

check whether all BAM read contain defined read groups? ( https://bioinformatics.stackexchange.com/questions/2590/ )

```
java -jar dist/samjdk.jar -e 'return record.getReadGroup()==null;'  input.bam
```

### Example 5

select chimeric reads virus/host

```
java -jar dist/samjdk.jar -e '
    final Set<String> virus_chrom = new HashSet<>(Arrays.asList("viral_seg1", "viral_seg2"));
    if(record.getReadUnmappedFlag()) return false;
    if(record.getReadPairedFlag() && !record.getMateUnmappedFlag())
        {
        if( virus_chrom.contains(record.getReferenceName()) !=
            virus_chrom.contains(record.getMateReferenceName())
            )
            {
            return true;
            }
        }
    for(final SAMRecord other: SAMUtils.getOtherCanonicalAlignments(record))
        {
        if( virus_chrom.contains(record.getReferenceName()) !=
            virus_chrom.contains(other.getReferenceName())
            )
            {
            return true;
            }
        }
    return false;' input.bam
```

### Example

cigar string with deletion >= 1kb

```
$ java -jar dist/samjdk.jar -e 'return !record.getReadUnmappedFlag() && record.getCigar().getCigarElements().stream().anyMatch(C->C.getLength()>=1000 && (C.getOperator()==CigarOperator.N || C.getOperator()==CigarOperator.D));'  in.bam
```

### Example

in *PAIRED* mode find some **pairs** of reads where at least one read startsWith `AAAAAAG`. See [https://bioinformatics.stackexchange.com/questions/2812](https://bioinformatics.stackexchange.com/questions/2812).


```
java -jar picard.jar SortSam I=S1.bam O=/dev/stdout  SO=queryname |\
java -jar dist/samjdk.jar --samoutputformat BAM --pair -e 'return records.stream().anyMatch(R->R.getReadString().startsWith("AAAAAAG"));' |\
samtools sort -T tmp - | samtools view -h - 
```

output:

```
@HD	VN:1.5	GO:none	SO:coordinate
@SQ	SN:rotavirus	LN:1074
@RG	ID:S1	SM:S1
@PG	ID:bwa	PN:bwa	VN:0.7.12-r1044	CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_01_R1.fq.gz S1_01_R2.fq.gz
@PG	ID:bwa.1	PN:bwa	VN:0.7.12-r1044	CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_02_R1.fq.gz S1_02_R2.fq.gz
@PG	ID:bwa.2	PN:bwa	VN:0.7.12-r1044	CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_03_R1.fq.gz S1_03_R2.fq.gz
rotavirus_34_627_6:0:0_4:0:0_b6	99	rotavirus	34	60	70M	=	558	594	GATGGATTCTCCTCAGCAGATGGTAAGGTCTATTATAAAACCTTCTTTTGAAGCTGCAGTTGTTGCTGCT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:6G3A16C8T2T0A29	RG:Z:S1	NM:i:6	AS:i:40	XS:i:0
rotavirus_40_627_9:0:0_3:0:0_ab	99	rotavirus	40	42	61M9S	=	558	588	GTCTACTCAGCAGATTGTAATCTCTCTTAATCATACTTCATTTGAAGCTGCCGTTGTTGCTTCTCCTTCA	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:15G4G4A3T1A7T11A9	RG:Z:S1	NM:i:7	AS:i:26	XS:i:19
rotavirus_52_627_5:0:0_5:0:0_2f8	163	rotavirus	52	60	70M	=	558	576	GATTGTAAGCTCTAATATTAAAACTTCTTTAGAAGGTGCAGTTGTTGCTGCTACTTCAACATTAGAATTA	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:3G10T6T8T4C34	RG:Z:S1	NM:i:5	AS:i:46	XS:i:0
rotavirus_132_627_6:0:0_6:0:0_1ee	163	rotavirus	132	60	70M	=	558	496	AATATGATTACAATGAAGTCTTTACCAGAGTTAAAAGTAAATATGCTTATGTGCTGTAAGACTCTGGTGT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:19A22T2A7A2G1T11	RG:Z:S1	NM:i:6	AS:i:40	XS:i:0
rotavirus_160_646_2:0:0_6:0:0_2e2	99	rotavirus	160	60	70M	=	577	487	AGTTAAAAGTAAATTTGATTCTGTGATGGATGACACTGGTGTTAAAAACAATCTTTTGGGTAAAGCTATA	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:20A13T35	RG:Z:S1	NM:i:2	AS:i:60	XS:i:0
rotavirus_239_683_6:0:0_4:0:0_a7	163	rotavirus	239	60	64M6S	=	614	445	CAGGCGTTAAATGGAAAGATTAGCTCAGCTATAAGAAATAGCAATTGGATGAGTGATTCTAAAAGGGTAG	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:18T13T8A10C11	RG:Z:S1	NM:i:4	AS:i:44	XS:i:0
rotavirus_132_627_6:0:0_6:0:0_1ee	83	rotavirus	558	60	70M	=	132	-496	AAAAAAGATTTGAATAACTATAACAGAGGGTTATTGACAAATACAATTCGTGGGTACAAAAAGCGAAGAA	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:15C4A12A3G9A1T20	RG:Z:S1	NM:i:6	AS:i:40	XS:i:0
rotavirus_34_627_6:0:0_4:0:0_b6	147	rotavirus	558	60	70M	=	34	-594	AAAAAAGATTTGAATCACTAAAACTGAGGGTTAATGAGAAATACAATACGTCGGTACATAAAGCGAAGAA	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:24A24T1G6A11	RG:Z:S1	NM:i:4	AS:i:50	XS:i:0
rotavirus_40_627_9:0:0_3:0:0_ab	147	rotavirus	558	60	70M	=	40	-588	AAAAAAGATTTGAATCACTAATACAGAGGGTTAATGAGATATACAATACTTGGGTACAAAAACCGAAGAA	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:21A17A22G7	RG:Z:S1	NM:i:3	AS:i:55	XS:i:0
rotavirus_52_627_5:0:0_5:0:0_2f8	83	rotavirus	558	60	70M	=	52	-576	AAAAAAGATTTGAATCACGATAACAGAGGGTGAATGAGAAATACATTACTTCGGTACAAAAAGCGAAGAA	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:18T1A10T13A5G18	RG:Z:S1	NM:i:5	AS:i:45	XS:i:0
rotavirus_160_646_2:0:0_6:0:0_2e2	147	rotavirus	577	60	70M	=	160	-487	AAAAAAGAGGGGTTATGAGTAATACAATACTTGGGTACAAAAAGAGATGAAAGTAAATGAAAATATGTAC	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:4C6T1A5A24C2A22	RG:Z:S1	NM:i:6	AS:i:41	XS:i:0
rotavirus_239_683_6:0:0_4:0:0_a7	83	rotavirus	614	60	70M	=	239	-445	AAAAAAGCGAAGAAAGTAAATGTAAATAGGTACTCTCTTCAGAATGTTATCTCACAACAGCAAAAACAAA	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:0C21A5T36C4	RG:Z:S1	NM:i:4	AS:i:54	XS:i:0
```

### Example

> Question: "How to extract reads with a known variant form a bam file"

Searching for a base 'T' on 'rotavirus' at 1044.

```java
final String contig= "rotavirus";
final int mutpos = 1044;
final char mutbase='T';
if(record.getReadUnmappedFlag()) return false;
if(!record.getContig().equals(contig)) return false;
if(record.getEnd() < mutpos) return false;
if(record.getStart() > mutpos) return false;
int readpos = record.getReadPositionAtReferencePosition(mutpos);
if(readpos<1) return false;
readpos--;
final byte[]    bases= record.getReadBases();
if(bases[readpos]==mutbase) return true;
return false;
```

### Example

Remove  Double clipped reads

```
java -jar dist/samjdk.jar -e 'if(record.getReadUnmappedFlag()) return true;final Cigar c=record.getCigar();if(c==null || c.numCigarElements()<2) return true; return !(c.getFirstCigarElement().getOperator().isClipping() && c.getLastCigarElement().getOperator().isClipping()) ;' input.bam
```
### Example

get discordant reads

```
$ java -jar dist/samjdk.jar -e 'return record.getReadPairedFlag() && !record.getReadUnmappedFlag() && !record.getMateUnmappedFlag() && !record.getReferenceName().equals(record.getMateReferenceName());'  in.bam
```


### Example

exclude clipped reads:

```
java -jar dist/samjdk.jar -e 'return record.getReadUnmappedFlag() || record.getCigar()==null || !record.getCigar().isClipped();'
```

###  Example

BAM file for reads that map 100% identity

```
java -jar  dist/samjdk.jar -e 'return !record.getReadUnmappedFlag() && record.getCigarString().equals("100M") &&  (record.getIntegerAttribute("NM")==null || record.getIntegerAttribute("NM").intValue()==0);' in.bam
```

### Example

> .. modify the vcf so that a snp that is detected where the read depth is more than the average read depth at the sample will change to the reference allele.

```
VariantContextBuilder vcb=new VariantContextBuilder(variant);
vcb.rmAttribute("AC");
vcb.rmAttribute("AF");
final double avgdp = variant.getGenotypes().stream().filter(G->G.hasDP()).mapToInt(G->G.getDP()).average().getAsDouble();
vcb.genotypes(variant.getGenotypes().
    stream().
    map(G->{
        if(!G.isCalled()) return G;
        if(!G.hasDP()) return G;
        if((double)G.getDP()> avgdp) return G;
        final List<Allele> aL=new ArrayList<>();
        while(aL.size()<G.getPloidy()) aL.add(variant.getReference());
        return new GenotypeBuilder(G).alleles(aL).make();
        }).
    collect(Collectors.toList())
    );
return vcb.make();
```

### Example:

https://bioinformatics.stackexchange.com/questions/3565/

> Subset smaller BAM to contain several thousand rows from multiple chromosomes

```
$ java -jar dist/samjdk.jar --body -e \
  'Map<String,Integer> c=new HashMap<>(); public Object apply(SAMRecord r) {int n=c.getOrDefault(r.getContig(),0);if(n>=5000) return false; c.put(r.getContig(),n+1); return r;}' \
 input.bam
```

### Example:

>  getting a list of alignment with a particular bam tag from a bam file

```bash
$ java -jar dist/samjdk.jar --body \
     -e 'Set<String> xc = null;  public Object apply( SAMRecord R) { if(xc==null) try {xc=new HashSet<>(IOUtil.slurpLines(new java.io.File("needed.txt")));} catch(Exception e) {throw new RuntimeIOException(e);} String att= R.getStringAttribute("XC"); return att!=null && xc.contains(att);}' \
     input.bam
```

### Example:

> remove reads where clipping > 33%

```bash
$ java -jar dist/samjdk.jar -e 'return record.getReadUnmappedFlag() || record.getCigar()==null || record.getCigar().getCigarElements().stream().filter(C->C.getOperator().isClipping()).mapToInt(C->C.getLength()).sum() / (double)record.getCigar().getReadLength() < 0.33;' input.bam
```

