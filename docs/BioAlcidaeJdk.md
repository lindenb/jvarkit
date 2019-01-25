# BioAlcidaeJdk

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

java-based version of awk for bioinformatics


## Usage

```
Usage: bioalcidaejdk [options] Files
  Options:
    --body
      user's code is the whole body of the filter class, not just the 'apply' 
      method. 
      Default: false
    -e, --expression
      inline java expression
    -F, --format
      force format: one of VCF BAM SAM FASTQ
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --import
      [20180312] add/import those java packages/classes in the code header. 
      Multiple separated by space/colon/comma. .eg: 'java.util.StringBuilder 
      java.awt.*' . 	Useful if those packages are not already defined in the 
      default code.
      Default: <empty string>
    --nocode
       Don't show the generated code
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -f, --scriptfile
      java body file
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * vcf
 * java
 * jdk



## See also in Biostars

 * [https://www.biostars.org/p/264894](https://www.biostars.org/p/264894)
 * [https://www.biostars.org/p/275714](https://www.biostars.org/p/275714)
 * [https://www.biostars.org/p/279535](https://www.biostars.org/p/279535)
 * [https://www.biostars.org/p/279942](https://www.biostars.org/p/279942)
 * [https://www.biostars.org/p/284852](https://www.biostars.org/p/284852)
 * [https://www.biostars.org/p/285803](https://www.biostars.org/p/285803)
 * [https://www.biostars.org/p/288324](https://www.biostars.org/p/288324)
 * [https://www.biostars.org/p/293237](https://www.biostars.org/p/293237)
 * [https://www.biostars.org/p/295040](https://www.biostars.org/p/295040)
 * [https://www.biostars.org/p/297983](https://www.biostars.org/p/297983)
 * [https://www.biostars.org/p/299255](https://www.biostars.org/p/299255)
 * [https://www.biostars.org/p/304780](https://www.biostars.org/p/304780)
 * [https://www.biostars.org/p/305174](https://www.biostars.org/p/305174)
 * [https://www.biostars.org/p/305743](https://www.biostars.org/p/305743)
 * [https://www.biostars.org/p/308310](https://www.biostars.org/p/308310)
 * [https://www.biostars.org/p/308554](https://www.biostars.org/p/308554)
 * [https://www.biostars.org/p/309013](https://www.biostars.org/p/309013)
 * [https://www.biostars.org/p/311363](https://www.biostars.org/p/311363)
 * [https://www.biostars.org/p/298361](https://www.biostars.org/p/298361)
 * [https://www.biostars.org/p/324900](https://www.biostars.org/p/324900)
 * [https://www.biostars.org/p/326294](https://www.biostars.org/p/326294)
 * [https://www.biostars.org/p/326765](https://www.biostars.org/p/326765)
 * [https://www.biostars.org/p/329423](https://www.biostars.org/p/329423)
 * [https://www.biostars.org/p/330752](https://www.biostars.org/p/330752)
 * [https://www.biostars.org/p/334253](https://www.biostars.org/p/334253)
 * [https://www.biostars.org/p/335056](https://www.biostars.org/p/335056)
 * [https://www.biostars.org/p/335692](https://www.biostars.org/p/335692)
 * [https://www.biostars.org/p/336206](https://www.biostars.org/p/336206)
 * [https://www.biostars.org/p/338031](https://www.biostars.org/p/338031)
 * [https://www.biostars.org/p/356474](https://www.biostars.org/p/356474)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bioalcidaejdk
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bioalcidae/BioAlcidaeJdk.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bioalcidae/BioAlcidaeJdk.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bioalcidae/BioAlcidaeJdkTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bioalcidae/BioAlcidaeJdkTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bioalcidaejdk** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

 * "bioalcidae, samjs and vcffilterjs: object-oriented formatters and filters for bioinformatics files" . Bioinformatics, 2017. Pierre Lindenbaum & Richard Redon  [https://doi.org/10.1093/bioinformatics/btx734](https://doi.org/10.1093/bioinformatics/btx734).


Bioinformatics file java-based reformatter. Something like awk for VCF, BAM, SAM...

This program takes as input a VCF or a BAM on stdin or as a file.
The user provides a piece of java code that will be compiled at runtime an executed.

## Why  this name, 'BioAlcidae' ?

As 'bioalcidae' looks like an 'awk' for bioinformatics, we used '[Alcidae](https://en.wikipedia.org/wiki/Alcidae)', the taxonomic Family of the '[auk](https://en.wikipedia.org/wiki/Auk)' species.

## History

  * 2019-01 migrating to openjdk11: switched to in-memory compiling to external compiling.

## Base classes 

At the time of writing, we have:

```java
    public static abstract class AbstractHandler
    	{
    	// output file
    	protected PrintStream out = System.out;
    	// input file name
    	protected String inputFile = null;
    	// called at begin
    	public void initialize() {}
    	// called at end
    	public void dispose() {}
    	// users MUST implement the body of that function
    	public abstract void execute() throws Exception;
    	
    	public void print(final Object o) { this.out.print(o);}
    	public void println() { this.out.println();}
    	public void println(final Object o) { this.print(o);this.println();}

    	}
    
    public static abstract class VcfHandler extends AbstractHandler
		{
    	protected VcfTools tools = null;
    	protected VCFHeader header = null;
    	protected VCFIterator iter = null;
		public Stream<VariantContext> stream()
			{
			return StreamSupport.stream(
					new IterableAdapter<VariantContext>(this.iter).spliterator(),
					false);
			}		}

    public static abstract class SAMHandler extends AbstractHandler
		{
		protected SamReader in=null;
		protected SAMFileHeader header=null;
		protected SAMRecordIterator iter=null;
		public Stream<SAMRecord> stream()
			{
			return StreamSupport.stream(
					new IterableAdapter<SAMRecord>(this.iter).spliterator(),
					false);
			}
		}
    public static abstract class FastqHandler extends AbstractHandler
		{
		protected FastqReader iter=null;
		public Stream<FastqRecord> stream()
			{
			return StreamSupport.stream(
					new IterableAdapter<FastqRecord>(this.iter).spliterator(),
					false);
			}
		}
		
    public static abstract class FastaHandler extends AbstractHandler
		{
		protected CloseableIterator<FastaSequence> iter=null;
		public Stream<FastaSequence> stream()
			{
			return StreamSupport.stream(
					new IterableAdapter<FastaSequence>(this.iter).spliterator(),
					false);
			}
		}

```

## VCF 

when reading a VCF, a new class extending `VcfHandler` will be compiled. The user's code will be inserted as:

```
     1  import java.util.*;
     2  import java.util.stream.*;
     3  import java.util.function.*;
     4  import htsjdk.samtools.*;
     5  import htsjdk.variant.variantcontext.*;
     6  import htsjdk.variant.vcf.*;
     7  import javax.annotation.processing.Generated;
     8  @Generated(value="BioAlcidaeJdk",date="2017-07-12T10:00:49+0200")
     9  public class Custom1694491176 extends com.github.lindenb.jvarkit.tools.bioalcidae.BioAlcidaeJdk.VcfHandler {
    10    public Custom1694491176() {
    11    }
    12    @Override
    13    public void execute() throws Exception {
    14    
	15     // user's code is inserted here <=======================
    16     
    17     }
    18  }
```

## SAM

when reading a SAM/BAM, a new class extending `SAMHandler` will be compiled. The user's code will be inserted as:


```java
     1  import java.util.*;
     2  import java.util.stream.*;
     3  import java.util.function.*;
     4  import htsjdk.samtools.*;
     5  import htsjdk.variant.variantcontext.*;
     6  import htsjdk.variant.vcf.*;
     7  import javax.annotation.processing.Generated;
     8  @Generated(value="BioAlcidaeJdk",date="2017-07-12T10:09:20+0200")
     9  public class Custom1694491176 extends com.github.lindenb.jvarkit.tools.bioalcidae.BioAlcidaeJdk.SAMHandler {
    10    public Custom1694491176() {
    11    }
    12    @Override
    13    public void execute() throws Exception {
    14    
	15     // user's code is inserted here <=======================
    16     
    17     }
    18  }
```


## FASTQ

when reading a Fastq, a new class extending `FastqHandler` will be compiled. The user's code will be inserted as:


```java
 1  import java.util.*;
 2  import java.util.stream.*;
 3  import java.util.function.*;
 4  import htsjdk.samtools.*;
 5  import htsjdk.samtools.util.*;
 6  import htsjdk.variant.variantcontext.*;
 7  import htsjdk.variant.vcf.*;
 8  import javax.annotation.processing.Generated;
 9  @Generated(value="BioAlcidaeJdk",date="2017-07-12T10:46:47+0200")
10  public class BioAlcidaeJdkCustom220769712 extends com.github.lindenb.jvarkit.tools.bioalcidae.BioAlcidaeJdk.FastqHandler {
11    public BioAlcidaeJdkCustom220769712() {
12    }
13    @Override
14    public void execute() throws Exception {
15    
16     // user's code is inserted here <=======================
17     
18     }
19  }

```

## Fasta

when reading a Fasta, a new class extending `FastaHandler` will be compiled. The user's code will be inserted as:

```
 1  import java.util.*;
 2  import java.util.stream.*;
 3  import java.util.function.*;
 4  import htsjdk.samtools.*;
 5  import htsjdk.samtools.util.*;
 6  import htsjdk.variant.variantcontext.*;
 7  import htsjdk.variant.vcf.*;
 8  import com.github.lindenb.jvarkit.util.bio.fasta.FastaSequence;
 9  import javax.annotation.processing.Generated;
10  @Generated(value="BioAlcidaeJdk",date="2017-07-12T14:26:39+0200")
11  public class BioAlcidaeJdkCustom298960668 extends com.github.lindenb.jvarkit.tools.bioalcidae.BioAlcidaeJdk.FastaHandler {
12    public BioAlcidaeJdkCustom298960668() {
13    }
14    @Override
15    public void execute() throws Exception {
16     // user's code starts here 
17     
18      //user's code ends here 
19     }
20  }

```



## Examples

### Example:

count the SNP having a ID in a VCF file:

```
$ java -jar dist/bioalcidaejdk.jar -e 'println(stream().filter(V->V.isSNP()).filter(V->V.hasID()).count());' input.vcf.gz

953
```

### Example:

count the mapped SAM Reads in a BAM file:

```
$ java -jar dist/bioalcidaejdk.jar -e 'println(stream().filter(R->!R.getReadUnmappedFlag()).count());' input.bam
5518 
```

### Example:

count the  Reads starting with "A' in a FASTQ file:

```
$ java -jar dist/bioalcidaejdk.jar -e 'println(stream().filter(R->R.getReadString().startsWith("A")).count());' in.fastq.gz
218 
```

## Example

generate a matrix: first row is sample names, one row per variant, genotypes (0=HOM_REF,1=HET,2=HOM_VAR,other=-9)

```
java -jar dist/bioalcidaejdk.jar -e 'println(String.join(",",header.getSampleNamesInOrder())); stream().forEach(V->println(V.getGenotypes().stream().map(G->{if(G.isHomRef()) return "0";if(G.isHet()) return "1"; if(G.isHomVar()) return "2"; return "-9";}).collect(Collectors.joining(","))));' input.vcf
```

## Example

longest fasta  length

```
/java -jar dist/bioalcidaejdk.jar -e 'println(stream().mapToInt(S->S.length()).max().getAsInt());' input.fasta
```

## Example

Which tool to calculate per site stats on vcf file?

```
java -jar dist/bioalcidaejdk.jar -e 'print("POS\t");for(GenotypeType GT : GenotypeType.values()) print("\t"+GT); println(); stream().forEach(V->{print(V.getContig()+":"+V.getStart()+":"+V.getReference().getDisplayString());for(GenotypeType GT : GenotypeType.values()) print("\t"+V.getGenotypes().stream().filter(G->G.getType()==GT).count()); println();});' in.vcf | column -t

POS               NO_CALL  HOM_REF  HET  HOM_VAR  UNAVAILABLE  MIXED
rotavirus:51:A    0        3        0    1        0            0
rotavirus:91:A    0        3        1    0        0            0
rotavirus:130:T   0        3        1    0        0            0
(...)
```

## Example

BED from cigar string with deletion >= 1kb

```
java -jar dist/bioalcidaejdk -e  'stream().filter(R->!R.getReadUnmappedFlag()).forEach(R->{ int refpos = R.getStart(); for(CigarElement ce:R.getCigar()) { CigarOperator op = ce.getOperator(); int len = ce.getLength(); if(len>=1000 && (op.equals(CigarOperator.N) || op.equals(CigarOperator.D))) { println(R.getContig()+"\t"+(refpos-1)+"\t"+(refpos+len)); } if(op.consumesReferenceBases())  refpos+=len; } }); ' in.bam
```

## Example

Printing the *aligned* reads in a given region of the reference (indels ignored)

```java
final String chrom = "rotavirus";
final int start=150;
final int end=200;

stream().filter(R->!R.getReadUnmappedFlag()).
    filter(R->R.getContig().equals(chrom) && !(R.getEnd()<start || R.getStart()>end)).
    forEach(R->{

    for(int pos=start; pos< end ;++pos)
        {
        char base='.';
        int ref=R.getStart();
        int read=0;
        for(CigarElement ce:R.getCigar())
            {
            CigarOperator op=ce.getOperator();
            switch(op)
                {
                case P: break;
                case H: break;
                case S : read+=ce.getLength(); break;
                case D: case N: ref+=ce.getLength();break;
                case I: read+=ce.getLength();break;
                case EQ:case X: case M:
                    {
                    for(int i=0;i< ce.getLength() ;i++)
                        {
                        if(ref==pos)
                            {
                            base = R.getReadString().charAt(read);
                            break;
                            }
                        if(ref>pos) break;
                        read++;
                        ref++;
                        }
                    break;
                    }
                default:break;
                }
            if(base!='.')break;
            }
        print(base);
        }
    println();
});
```

## Example

print information about reads overlaping the base 'rotavirus:200' and the tag 'XS'


```
$ samtools view -h ~/src/gatk-ui/testdata/S1.bam "rotavirus:200-200" | java -jar dist/bioalcidaejdk.jar -F SAM -e 'stream().forEach(record->{final String contig= "rotavirus"; final int mutpos = 200; if(record.getReadUnmappedFlag()) return ;if(!record.getContig().equals(contig)) return ;if(record.getEnd() < mutpos) return ;if(record.getStart() > mutpos) return ; int readpos = record.getReadPositionAtReferencePosition(mutpos); if(readpos<1) return ; readpos--; final byte[] bases= record.getReadBases();println(record.getReadName()+" "+record.getContig()+" "+record.getStart()+" "+record.getEnd()+" base["+mutpos+"]="+(char)bases[readpos]+" XS="+record.getAttribute("XS"));});'
```

## Example


```
$ java -jar dist/bioalcidaejdk.jar -e 'stream().forEach(V->{println(V.getContig()+"\t"+V.getStart()+"\t"+V.getReference().getDisplayString()+"\t"+V.getAlternateAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(","))+"\t"+V.getGenotypes().stream().map(G->G.getSampleName()+"="+G.getAlleles().stream().filter(A->A.isCalled()).map(A->A.getDisplayString()).collect(Collectors.joining("/"))).collect(Collectors.joining("\t")));});' src/test/resources/test_vcf01.vcf

1   956852  C   T   S1= S2=T/T  S3=C/T  S4=C/T  S5=T/T  S6=C/T
1   959155  G   A   S1= S2=G/A  S3=G/A  S4=G/A  S5=A/A  S6=G/A
1   959169  G   C   S1= S2=G/C  S3=G/C  S4=G/C  S5=C/C  S6=G/C
1   959231  G   A   S1= S2=G/A  S3=G/A  S4=G/A  S5=A/A  S6=G/A
1   960409  G   C   S1= S2=G/C  S3=G/C  S4=G/C  S5=C/C  S6=G/C
1   962210  A   G   S1=A/G  S2=A/G  S3=A/G  S4=A/G  S5= S6=A/G
```

## Example

'Merging of sequence chunks' [https://www.biostars.org/p/288324/#288466](https://www.biostars.org/p/288324/#288466)

```java
final Map<String,List<FastaSequence>> id2seqs= new HashMap<>();
stream().forEach(S->{
    String n = S.getName();
    if(!id2seqs.containsKey(n)) id2seqs.put(n,new ArrayList<>());
    id2seqs.get(n).add(S);
    });

for(final List<FastaSequence> seqs: id2seqs.values())
{
println(">"+seqs.get(0).getName());
int len = seqs.stream().mapToInt(S->S.length()).max().getAsInt();
for(int i=0;i< len;i++)
    {
    char c='\0';
    for(FastaSequence s:seqs)
        {
        if(i>=s.length() || s.charAt(i)=='-'  ||  c==s.charAt(i)) continue;
        if(c=='\0')
            {
            c=s.charAt(i);
            }
        else
            {   
            c='X';
            }
        }
    print(c=='\0'?'-':c);
    }
println();
}
```
## Example

> SNP,INDEL counting per chromosomes in vcf

```java
stream().
  map(V->V.getType().name()+" "+V.getContig()).
  collect(Collectors.groupingBy(Function.identity(), Collectors.counting())).forEach((K,V)->{println(K+" : "+V);});
```

## Example

> Print Minimum Depth from a VCF

```
$ java -jar dist/bioalcidaejdk.jar -e 'stream().forEach(V->{println(V.getContig()+":"+V.getStart()+":"+V.getReference().getDisplayString()+"\t"+V.getGenotypes().stream().filter(G->G.hasDP()).mapToInt(G->G.getDP()).min().orElse(-1));});' input.vcf
```

## Example


>  remove sequences with non-canonical nucleotides from fasta file

```
$ java -jar dist/bioalcidaejdk.jar -e 'stream().filter(F->java.util.regex.Pattern.matches("^[ATGCatgc]+$",F)).forEach(S->println(">"+S.getName()+"\n"+S));' input.fa
```

## Example

> Parse BAM by insertion size and get genomic coordinates

https://bioinformatics.stackexchange.com/questions/3492/parse-bam-by-insertion-size-and-get-genomic-coordinates

```
$ java -jar dist/bioalcidaejdk.jar -e 'stream().filter(R->!R.getReadUnmappedFlag() && R.getCigar().getCigarElements().stream().anyMatch(C->C.getLength()>100 && (C.getOperator().equals(CigarOperator.N) || C.getOperator().equals(CigarOperator.D)))).forEach(R->println(R.getContig()+"\t"+R.getStart()+"\t"+R.getEnd()));' input.bam
```

## Example

Print BED file with NM attributes

```
java -jar dist/bioalcidaejdk.jar -e 'stream().filter(R->!R.getReadUnmappedFlag() && R.hasAttribute("NM")).forEach(R->println(R.getContig()+"\t"+(R.getStart()-1)+"\t"+R.getEnd()+"\t"+R.getIntegerAttribute("NM")));'
```

## Example

 Extract all genotype counts from phased data in vcf files

```
$ wget -q -O - "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" |\
gunzip -c |\
java -jar dist/bioalcidaejdk.jar -F VCF -e 'stream().forEach(V->println(V.getContig()+" "+V.getStart()+" " +V.getReference().getDisplayString()+" "+V.getGenotypes().stream().map(G->G.getAlleles().stream().map(A->String.valueOf(V.getAlleleIndex(A))).collect(Collectors.joining(G.isPhased()?"|":"/"))).collect(Collectors.groupingBy(Function.identity(), Collectors.counting()))));'

22 16050075 A {0|0=2503, 0|1=1}
22 16050115 G {0|0=2472, 1|0=16, 0|1=16}
22 16050213 C {0|0=2467, 1|0=18, 0|1=18, 1|1=1}
22 16050319 C {0|0=2503, 1|0=1}
22 16050527 C {0|0=2503, 0|1=1}
22 16050568 C {0|0=2502, 0|1=2}
22 16050607 G {0|0=2499, 1|0=3, 0|1=2}
22 16050627 G {0|0=2502, 1|0=2}
22 16050646 G {0|0=2503, 1|0=1}
22 16050654 A {0|0=1834, 1|0=3, 0|1=6, 0|2=48, 2|0=37, 0|3=255, 3|0=258, 4|0=6, 0|4=12, 2|3=1, 3|2=1, 3|3=41, 4|3=2}
22 16050655 G {0|0=2503, 1|0=1}
22 16050678 C {0|0=2502, 0|1=2}
22 16050679 G {0|0=2503, 1|0=1}
22 16050688 C {0|0=2503, 1|0=1}
22 16050732 C {0|0=2503, 0|1=1}
22 16050739 TA {0|0=2467, 1|0=18, 0|1=18, 1|1=1}
22 16050758 T {0|0=2503, 0|1=1}
22 16050783 A {0|0=2466, 0|1=13, 1|0=24, 1|1=1}
22 16050840 C {0|0=2478, 1|0=11, 0|1=15}
22 16050847 T {0|0=2498, 1|0=2, 0|1=4}
```

## Example

> How to get percent identity from bam file ?

using the edit distance 'NM' 

```
$ java -jar dist/bioalcidaejdk.jar -e 'stream().map(R->R.getReadUnmappedFlag()?0:(int)(100.0*(R.getReadLength()-R.getIntegerAttribute("NM"))/(double)R.getReadLength())).collect(Collectors.groupingBy(Function.identity(), Collectors.counting())).forEach((K,V)->{println(K+"\t"+V);});' input.bam
```


## Example

sliding windows of 1000 variant, shift by 500 variants

```
java -jar bioalcidaejdk.jar -e 'String prevContig=null; Consumer<List<VariantContext>> dump = (L)->{if(L.isEmpty()) return;System.out.println(L.get(0).getContig()+"\t"+(L.get(0).getStart()-1)+"\t"+L.get(L.size()-1).getEnd());}; final List<VariantContext> buffer=new ArrayList<>();while(iter.hasNext()) {VariantContext vc=iter.next();if(!vc.getContig().equals(prevContig)) {dump.accept(buffer);buffer.clear();prevContig=vc.getContig();} buffer.add(vc);if(buffer.size()>=1000) {dump.accept(buffer);buffer.subList(0,500).clear();}} dump.accept(buffer);' in.vcf
```

## number of homozygous and heterozygous non-reference from a multi sample VCF file?

```
java -jar dist/bioalcidaejdk.jar -e 'stream().flatMap(V->V.getGenotypes().stream()).filter(G->!G.isHomRef()).map(G->G.getSampleName()+"\t"+G.getType().name()).collect(Collectors.groupingBy(Function.identity(), Collectors.counting())).forEach((K,V)->println(K+"\t"+V));' src/test/resources/rotavirus_rf.vcf.gz
```

## Counting soft clipped bases and reads

```java
final long counts[]=new long[4];
stream().forEach(R->{
	counts[0]++;
	if(R.getReadUnmappedFlag()) return;
	Cigar cigar = R.getCigar();
	if(cigar==null || cigar.isEmpty()) return;
	counts[1]++;
	counts[2] += cigar.getReadLength();
	counts[3] += cigar.
		getCigarElements().
		stream().
		filter(C->C.getOperator().isClipping()).
		mapToInt(C->C.getLength()).
		sum();
	});

println("NUM-READS:"+counts[0]);
println("NUM-MAPPED-READS:"+counts[1]);
println("SUM-MAPPED-READS-LENGTH:"+counts[2]);
println("SUM-CLIPPING:"+counts[3]);
```



