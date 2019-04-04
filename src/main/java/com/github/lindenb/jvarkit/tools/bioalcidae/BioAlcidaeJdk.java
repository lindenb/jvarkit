/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
package com.github.lindenb.jvarkit.tools.bioalcidae;


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.lang.reflect.Constructor;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.OpenJdkCompiler;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.fasta.FastaSequence;
import com.github.lindenb.jvarkit.util.bio.fasta.FastaSequenceReader;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfTools;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Iso8601Date;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

/**

BEGIN_DOC

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



END_DOC
*/


@Program(name="bioalcidaejdk",
	description="java-based version of awk for bioinformatics",
	keywords={"sam","bam","vcf","java","jdk"},
	biostars={264894,275714,279535,279942,284852,285803,288324,293237,295040,
			297983,299255,304780,305174,305743,308310,308554,309013,311363,
			298361,324900,326294,326765,329423,330752,334253,335056,335692,336206,
			338031,356474},
	references="\"bioalcidae, samjs and vcffilterjs: object-oriented formatters and filters for bioinformatics files\" . Bioinformatics, 2017. Pierre Lindenbaum & Richard Redon  [https://doi.org/10.1093/bioinformatics/btx734](https://doi.org/10.1093/bioinformatics/btx734)."
	)
public class BioAlcidaeJdk
	extends Launcher
	{
	private static final Logger LOG = Logger.build(BioAlcidaeJdk.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-F","--format"},description="force format: one of VCF BAM SAM FASTQ")
	private String formatString = null;

    @Parameter(names={"-f","--scriptfile"},description="java body file")
    private File scriptFile=null;
    @Parameter(names={"-e","--expression"},description="inline java expression")
    private String scriptExpr=null;
	@Parameter(names={"--nocode"},description=" Don't show the generated code")
	private boolean hideGeneratedCode=false;
	@Parameter(names={"--body"},description="user's code is the whole body of the filter class, not just the 'apply' method.")
	private boolean user_code_is_body=false;
	@Parameter(names={"--import"},description="[20180312] add/import those java packages/classes in the code header. "
			+ "Multiple separated by space/colon/comma. .eg: 'java.util.StringBuilder java.awt.*' . "
			+ "	Useful if those packages are not already defined in the default code.")
	private String extraImport = "";

	
	@SuppressWarnings("unused")
	private static final Counter<?> _fool_javac = null;
	
    public static abstract class AbstractHandler
    	{
    	protected PrintStream out = System.out;
    	protected String inputFile = null;
    	public void initialize() {}
    	public void dispose() {}
    	public void print(final Object o) { this.out.print(o);}
    	public void println() { this.out.println();}
    	public void println(final Object o) { this.print(o);this.println();}
    	public abstract void execute() throws Exception;
    	}
    
    public static abstract class AbstractHandlerFactory<H extends AbstractHandler>
    	{
    	private File scriptFile  = null;
    	private String scriptExpr = null ;
    	private boolean user_code_is_body=false;
    	private boolean hideGeneratedCode = false;
    	private Constructor<H> ctor=null;
    	private Set<String> extraImportSet = new HashSet<>();

    	
    	public abstract int execute(final String inputFile,final PrintStream out) throws Exception;
    	protected abstract Class<H> getHandlerClass();
    	
    	
    	protected  BufferedReader openBufferedReader(final String inOrNull) throws IOException {
    		return(inOrNull==null?
    				new BufferedReader(new InputStreamReader(System.in)):
    				IOUtils.openURIForBufferedReading(inOrNull)
    				);
    		}
    	
    	@SuppressWarnings("unchecked")
		public Constructor<H> getConstructor() {
    		if(this.ctor !=null) return this.ctor;
    		if(this.scriptFile!=null && !StringUtil.isBlank(this.scriptExpr))
				{
				throw new JvarkitException.UserError("script file and expression both defined");
				}
			
			if(this.scriptFile==null && StringUtil.isBlank(this.scriptExpr))
				{
				throw new JvarkitException.UserError("script file or expression missing");
				}
			
			try {
				
				final Random rand= new  Random(System.currentTimeMillis());
				final String javaClassName =BioAlcidaeJdk.class.getSimpleName()+
						"Custom"+ Math.abs(rand.nextInt());
				
				final String baseClass = getHandlerClass().getName().replace('$', '.');
				final String code;
				
				if(this.scriptFile!=null)
					{
					code = IOUtil.slurp(this.scriptFile);
					}
				else
					{
					code = this.scriptExpr;
					}
				
				final String generatedClass = OpenJdkCompiler.getGeneratedAnnotationClassName();
				final StringWriter codeWriter=new StringWriter();
				final PrintWriter pw = new PrintWriter(codeWriter);
				pw.println("import java.util.*;");
				pw.println("import java.util.stream.*;");
				pw.println("import java.util.regex.*;");
				pw.println("import java.util.function.*;");
				pw.println("import htsjdk.samtools.*;");
				pw.println("import htsjdk.samtools.util.*;");
				pw.println("import htsjdk.variant.variantcontext.*;");
				pw.println("import htsjdk.variant.vcf.*;");
				pw.println("import com.github.lindenb.jvarkit.util.bio.fasta.FastaSequence;");
				pw.println("import htsjdk.variant.vcf.*;");
				
				pw.println("/** begin user's packages */");
				for(final String p:this.extraImportSet)
					{
					if(StringUtil.isBlank(p)) continue;
					pw.print("import ");
					pw.print(p);
					
					pw.println(";");
					
					}
				pw.println("/** end user's packages */");
				
				
				if(!StringUtil.isBlank(generatedClass)) {
					pw.println("@"+generatedClass+"(value=\""+BioAlcidaeJdk.class.getSimpleName()+"\",date=\""+ new Iso8601Date(new Date()) +"\")");
					}
				pw.println("public class "+javaClassName+" extends "+ baseClass +" {");
				
				pw.println("  public "+javaClassName+"() {");
				pw.println("  }");
				
				
				if(this.user_code_is_body)
					{
					pw.println("   // user's code starts here ");
					pw.println("   "+ code);
					pw.println("   // user's code ends here ");
					}
				else
					{
					pw.println("  @Override");
					pw.println("  public void execute() throws Exception {");
					pw.println("   // user's code starts here ");
					pw.println("   "+ code);
					pw.println("    //user's code ends here ");
					pw.println("   }");
					}
				pw.println("}");
				pw.flush();
				
				
				if(!this.hideGeneratedCode)
					{
					LOG.debug(" Compiling :\n" + OpenJdkCompiler.beautifyCode(codeWriter.toString()));
					}

				final OpenJdkCompiler inMemoryCompiler = OpenJdkCompiler.getInstance();
				final Class<?> compiledClass = inMemoryCompiler.compileClass(
						javaClassName,
						codeWriter.toString()
						);
				this.ctor = (Constructor<H>)compiledClass.getDeclaredConstructor();
				return this.ctor;
				}
			catch(final Throwable err) {
				throw new RuntimeException(err);
				}
			}
    	
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
			}		
		}

    
    public static class VcfHandlerFactory extends AbstractHandlerFactory<VcfHandler>
    	{
    	@Override
    	protected Class<VcfHandler> getHandlerClass() {
    		return VcfHandler.class;
    		}
    	@Override
    	public int execute(final String inputFile,final PrintStream out) throws Exception {
    		VcfHandler vcfHandler = null;
    		try 
	    		{
	    		vcfHandler= (VcfHandler)getConstructor().newInstance();
	    		vcfHandler.out = out;
	    		vcfHandler.inputFile = inputFile;
				//
				vcfHandler.iter = VCFUtils.createVCFIterator(inputFile);
				vcfHandler.header = vcfHandler.iter.getHeader();
				vcfHandler.tools = new VcfTools(vcfHandler.header);
				vcfHandler.initialize();
				vcfHandler.execute();
				return 0;
	    		}
    		catch(Throwable err)
    			{
    			LOG.error(err);
    			return -1;
    			}
    		finally
    			{
    			if(vcfHandler!=null) {
    				vcfHandler.dispose();
    				if(vcfHandler.out!=null) vcfHandler.out.flush();
    				CloserUtil.close(vcfHandler.out);
    				CloserUtil.close(vcfHandler.iter);
    				}
    			}
			}
    	}
    
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

    public static class SAMHandlerFactory extends AbstractHandlerFactory<SAMHandler>
		{
    	@Override
    	protected Class<SAMHandler> getHandlerClass() {
    		return SAMHandler.class;
    		}
		@Override
		public int execute(final String inputFile,final PrintStream out) throws Exception {
			SAMHandler samHandler= null;
			try
				{
				samHandler = 	(SAMHandler)this.getConstructor().newInstance();
				samHandler.out = out;
				samHandler.inputFile = inputFile;
				//
				final htsjdk.samtools.SamReaderFactory srf= htsjdk.samtools.SamReaderFactory.makeDefault().validationStringency(htsjdk.samtools.ValidationStringency.LENIENT);
				if(inputFile==null)
					{
					samHandler.in =  srf.open(htsjdk.samtools.SamInputResource.of(System.in));
					}
				else
					{
					samHandler.in = srf.open(htsjdk.samtools.SamInputResource.of(inputFile));
					}			
				samHandler.header = samHandler.in.getFileHeader();
				samHandler.iter = samHandler.in.iterator();	
				
				samHandler.initialize();
				samHandler.execute();
				return 0;
				}
			catch (final Throwable err) {
				LOG.error(err);
				return -1;
				}
			finally
				{
				if(samHandler!=null) 
					{
					samHandler.dispose();
					CloserUtil.close(samHandler.iter);
					CloserUtil.close(samHandler.in);
					if(samHandler.out!=null) samHandler.out.flush();
					CloserUtil.close(samHandler.out);
					samHandler=null;
					}
				}
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
    
    
    public static class FastqHandlerFactory extends AbstractHandlerFactory<FastqHandler>
		{
    	@Override
    	protected Class<FastqHandler> getHandlerClass() {
    		return FastqHandler.class;
    		}
		@Override
		public int execute(final String inputFile,final PrintStream out) throws Exception {
			FastqHandler fqHandler= null;
			try {
				fqHandler = (FastqHandler)this.getConstructor().newInstance();
				fqHandler.out = out;
				fqHandler.inputFile = inputFile;
				fqHandler.iter = new FastqReader(super.openBufferedReader(inputFile));
				
				fqHandler.initialize();
				fqHandler.execute();
				return 0;
				}
			catch (final Throwable err) {
				LOG.error(err);
				return -1;
				}
			finally
				{
				if(fqHandler!=null) 
					{
					fqHandler.dispose();
					CloserUtil.close(fqHandler.iter);
					if(fqHandler.out!=null) fqHandler.out.flush();
					CloserUtil.close(fqHandler.out);
					fqHandler=null;
					}
				}
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
    
    public static class FastaHandlerFactory extends AbstractHandlerFactory<FastaHandler>
		{
    	@Override
    	protected Class<FastaHandler> getHandlerClass() {
    		return FastaHandler.class;
    		}
    	@Override
    	public int execute(final String inputFile, final PrintStream out) throws Exception {
    		FastaHandler faHandler= null;
			try {
				faHandler = (FastaHandler)this.getConstructor().newInstance();
				faHandler.out = out;
				faHandler.inputFile = inputFile;
				//
				faHandler.iter =  new FastaSequenceReader().iterator(super.openBufferedReader(inputFile));
				
				faHandler.initialize();
				faHandler.execute();
				return 0;
				}
			catch (final Throwable err) {
				LOG.error(err);
				return -1;
				}
			finally
				{
				if(faHandler!=null) 
					{
					faHandler.dispose();
					CloserUtil.close(faHandler.iter);
					if(faHandler.out!=null) faHandler.out.flush();
					CloserUtil.close(faHandler.out);
					faHandler=null;
					}
				}
			}
		
		}

    public static abstract class SimpleLineHandler extends AbstractHandler
		{
		protected CloseableIterator<String> iter=null;
		public Stream<String> stream()
			{
			return StreamSupport.stream(
					new IterableAdapter<String>(this.iter).spliterator(),
					false);
			}
		}
    
    public static class SimpleLineHandlerHandlerFactory extends AbstractHandlerFactory<SimpleLineHandler>
		{
    	private static class LineIterator extends AbstractIterator<String> implements CloseableIterator<String> {
    		private final BufferedReader br;
    		LineIterator(final BufferedReader br) {
    			this.br = br;
    			}
    		@Override
    		protected String advance() {
    			try {
    				return br.readLine();
    				}
    			catch(Exception err) {
    				return null;
    			}
    			}
    		@Override
    		public void close() {
    			CloserUtil.close(br);
    			}
    	}
    	@Override
    	protected Class<SimpleLineHandler> getHandlerClass() {
    		return SimpleLineHandler.class;
    		}
    	@Override
    	public int execute(final String inputFile, final PrintStream out) throws Exception {
    		SimpleLineHandler lineHandler = null;
			try {
				lineHandler = (SimpleLineHandler)this.getConstructor().newInstance();
				lineHandler.out = out;
				lineHandler.inputFile = inputFile;
				//
				lineHandler.iter =  new LineIterator(IOUtils.openURIForBufferedReading(inputFile));
				
				lineHandler.initialize();
				lineHandler.execute();
				return 0;
				}
			catch (final Throwable err) {
				LOG.error(err);
				return -1;
				}
			finally
				{
				if(lineHandler!=null) 
					{
					lineHandler.dispose();
					CloserUtil.close(lineHandler.iter);
					if(lineHandler.out!=null) lineHandler.out.flush();
					CloserUtil.close(lineHandler.out);
					lineHandler=null;
					}
				}
			}
		
		}

    
	private enum FORMAT {
		VCF{
			@Override
			boolean canAs(final String src) {
				return src!=null && (Arrays.asList(IOUtil.VCF_EXTENSIONS).stream().anyMatch(EXT->src.endsWith(EXT)) );
			}
			},
		SAM{
			@Override
			boolean canAs(final String src) {
				return src!=null && (src.endsWith(".sam") || src.endsWith(".bam") );
			}
			},
		BAM{
			@Override
			boolean canAs(final String src) {
				return src!=null && (src.endsWith(".sam") || src.endsWith(".bam") );
			}},
		FASTQ{
			@Override
			boolean canAs(final String src) {
				return src!=null && (src.endsWith(".fq") || src.endsWith(".fastq") || src.endsWith(".fq.gz") || src.endsWith(".fastq.gz")  );
			}},
		FASTA{
			@Override
			boolean canAs(final String src) {
				return src!=null && (src.endsWith(".fa") || src.endsWith(".fasta") || src.endsWith(".fa.gz") || src.endsWith(".fasta.gz")  );
			}},
		TEXT {
			@Override
			boolean canAs(final String src) {
				return src!=null && (src.endsWith(".txt") || src.endsWith(".csv") || src.endsWith(".tsv") );
			}}
			;

		abstract boolean canAs(String src);
		};
		
		private FORMAT format= null;

	
	@Override
	public int doWork(final List<String> args) {
		AbstractHandlerFactory<?> abstractFactory = null;
		if(this.formatString!=null)
			{
			try {
				this.format=FORMAT.valueOf(this.formatString.toUpperCase());
				} catch (Exception err) {
				LOG.error(err);
				return -1;
				}
			}
		try
			{
			final String inputFile = oneFileOrNull(args);
			if(inputFile==null && this.format==null)
				{
				LOG.error("Format must be specified when input is stdin");
				return -1;
				}
			if(inputFile!=null && this.format==null)
				{
				for(final FORMAT t:FORMAT.values())
					{
					if(t.canAs(inputFile))
						{
						this.format=t;
						break;
						}
					}
				if(this.format==null)
					{
					LOG.error("Cannot get file format for "+inputFile+". Please specifiy using option -F");
					return -1;
					}
				}
			
			 
			
			switch(this.format)
				{
				case BAM: case SAM: abstractFactory = new SAMHandlerFactory(); break;
				case VCF: abstractFactory = new VcfHandlerFactory(); break;
				case FASTQ: abstractFactory = new FastqHandlerFactory(); break;
				case FASTA: abstractFactory = new FastaHandlerFactory(); break;
				case TEXT: abstractFactory = new SimpleLineHandlerHandlerFactory();break;
				default: throw new IllegalStateException("Not implemented: "+this.format);
				}
			abstractFactory.scriptExpr = this.scriptExpr;
			abstractFactory.scriptFile = this.scriptFile;
			abstractFactory.user_code_is_body = this.user_code_is_body ;
			abstractFactory.hideGeneratedCode = this.hideGeneratedCode ;
			abstractFactory.extraImportSet.addAll(
					Arrays.asList(this.extraImport.split("[ ,;]+"))
					);//blank will be ignored
			
			try
				{
				return abstractFactory.execute(inputFile,super.openFileOrStdoutAsPrintStream(this.outputFile));
				}
			catch(final Throwable err)
				{
				LOG.error(err);
				return -1;
				}
			finally
				{
				abstractFactory.ctor=null;
				}			
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			
			}
		}
	
	
	public static void main(String[] args) {
		new BioAlcidaeJdk().instanceMainWithExit(args);
	}
}
