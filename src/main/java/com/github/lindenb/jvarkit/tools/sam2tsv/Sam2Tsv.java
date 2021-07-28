/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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


History:
* 2014-11 : handle clipped bases
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.sam2tsv;

import java.io.PrintWriter;
import java.nio.file.Path;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.MultiBamLauncher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
/**

BEGIN_DOC


### Example
 


```
$  java -jar dist/sam2tsv.jar -R src/test/resources/toy.fa src/test/resources/toy.bam 


#Read-Name  Flag  MAPQ  CHROM  READ-POS0  READ-BASE  READ-QUAL  REF-POS1  REF-BASE  CIGAR-OP
r001        163   30    ref    0          T          .          7         T         M
r001        163   30    ref    1          T          .          8         T         M
r001        163   30    ref    2          A          .          9         A         M
r001        163   30    ref    3          G          .          10        G         M
r001        163   30    ref    4          A          .          11        A         M
r001        163   30    ref    5          T          .          12        T         M
r001        163   30    ref    6          A          .          13        A         M
r001        163   30    ref    7          A          .          14        A         M
r001        163   30    ref    8          A          .          .         .         I
r001        163   30    ref    9          G          .          .         .         I
r001        163   30    ref    10         A          .          .         .         I
r001        163   30    ref    11         G          .          .         .         I
r001        163   30    ref    12         G          .          15        G         M
r001        163   30    ref    13         A          .          16        A         M
r001        163   30    ref    14         T          .          17        T         M
r001        163   30    ref    15         A          .          18        A         M
r001        163   30    ref    .          .          .          19        G         D
r001        163   30    ref    16         C          .          20        C         M
r001        163   30    ref    17         T          .          21        T         M
r001        163   30    ref    18         G          .          22        G         M
r002        0     30    ref    0          A          .          8         T         S
r002        0     30    ref    1          A          .          .         .         I
r002        0     30    ref    2          A          .          .         .         I
r002        0     30    ref    3          A          .          9         A         M
r002        0     30    ref    4          G          .          10        G         M
r002        0     30    ref    5          A          .          11        A         M
r002        0     30    ref    6          T          .          12        T         M
r002        0     30    ref    7          A          .          13        A         M
r002        0     30    ref    8          A          .          14        A         M
r002        0     30    ref    .          .          .          .         .         P
r002        0     30    ref    9          G          .          .         .         I
r002        0     30    ref    .          .          .          .         .         P
r002        0     30    ref    10         G          .          .         .         I
r002        0     30    ref    11         G          .          15        G         M
r002        0     30    ref    12         A          .          16        A         M
(...)

```


## Example 2

sam2tsv can read data from a linux pipe.

```
samtools view -h input.bam | java -jar dist/sam2tsv.jar
```



### Citations


Sam2tsv was cited in : 

  * "Illumina TruSeq Synthetic Long-Reads Empower De Novo Assembly and Resolve Complex, Highly-Repetitive Transposable Elements" . McCoy RC, Taylor RW, Blauwkamp TA, Kelley JL, Kertesz M, et al. (2014) Illumina TruSeq Synthetic Long-Reads Empower De Novo Assembly and Resolve Complex, Highly-Repetitive Transposable Elements. PLoS ONE 9(9): e106689. doi: 10.1371/journal.pone.0106689  http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0106689
  * "High-Throughput Identification of Genetic Variation Impact on pre-mRNA Splicing Efficiency". Scott I Adamson, Lijun Zhan, Brenton R Graveley. doi: [https://doi.org/10.1101/191122](https://doi.org/10.1101/191122).
  * "Linkage of A-to-I RNA editing in metazoans and the impact on genome evolution "  Molecular Biology and Evolution, msx274, https://doi.org/10.1093/molbev/msx274
  * "Vex-seq: high-throughput identification of the impact of genetic variation on pre-mRNA splicing efficiency" Genome Biology201819:71 https://doi.org/10.1186/s13059-018-1437-x
  * "Accurate detection of m6A RNA modifications in native RNA sequences" Huanle Liu, Oguzhan Begik, Morghan C Lucas, Christopher E Mason, Schraga Schwartz, John S Mattick, Martin A Smith, Eva Maria Novoa bioRxiv 525741; doi: https://doi.org/10.1101/525741 
  * "DART-seq: an antibody-free method for global m6A detection"  Nature Methods https://doi.org/10.1038/s41592-019-0570-0
  * "Thiouridine-to-Cytidine Conversion Sequencing (TUC-Seq) to Measure mRNA Transcription and Degradation Rates" The Eukaryotic RNA Exosome. Nov 2019. https://doi.org/10.1007/978-1-4939-9822-7_10
  * "Evolutionary forces on A-to-I RNA editing revealed by sequencing individual honeybee drones". Yuange Duan, Shengqian Dou, Jiaxing Huang, Eli Eisenberg, Jian Lu . 2020 . https://doi.org/10.1101/2020.01.15.907287
  * "Sci-fate characterizes the dynamics of gene expression in single cells" (2020) Nat Biotechnol (2020). https://doi.org/10.1038/s41587-020-0480-9
  * "Mutations in virus-derived small RNAs" Nigam Deepti; LaTourrette, Katherine; Garcia-Ruiz, Hernan. Scientific Reports (Nature Publisher Group); London Vol. 10, Iss. 1,  (2020). DOI:10.1038/s41598-020-66374-2 
  * Qiu, Q., Hu, P., Qiu, X. et al. Massively parallel and time-resolved RNA sequencing in single cells with scNT-seq. Nat Methods (2020). https://doi.org/10.1038/s41592-020-0935-4
  * FUJIKURA, K. et al. Multiregion whole-exome sequencing of intraductal papillary mucinous neoplasms reveals frequent somatic KLF4 mutations predominantly in low-grade regions. Gut, [s. l.], 2020. DOI 10.1136/gutjnl-2020-321217
  * Gao, Y., Liu, X., Wu, B. et al. Quantitative profiling of N6-methyladenosine at single-base resolution in stem-differentiating xylem of Populus trichocarpa using Nanopore direct RNA sequencing. Genome Biol 22, 22 (2021). https://doi.org/10.1186/s13059-020-02241-7
  * Liu H., Begik O., Novoa E.M. (2021) EpiNano: Detection of m6A RNA Modifications Using Oxford Nanopore Direct RNA Sequencing. In: McMahon M. (eds) RNA Modifications. Methods in Molecular Biology, vol 2298. Humana, New York, NY. https://doi.org/10.1007/978-1-0716-1374-0_3

END_DOC
*/
@Program(name="sam2tsv",
	description="Prints the SAM alignments as a TAB delimited file.",
	keywords={"sam","bam","table","cram","tsv"},
	biostars={157232,59647,253828,264875,277493},
	modificationDate="20210304",
	creationDate="20170712"
	)
public class Sam2Tsv
	extends MultiBamLauncher
	{
	private static final Logger LOG = Logger.build(Sam2Tsv.class).make();

	@Parameter(names={"-o","--output"},description= OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-N","--skip-N"},description="Skip 'N' operator")
	private boolean skip_N=false;


	private ReferenceSequenceFile indexedFastaSequenceFile=null;
	private GenomicSequence genomicSequence=null;
	private PrintWriter out = null;
	
	// current row members
	private SAMRecord rec;
	private byte readbases[];
	private byte readQuals[];
	private int readPos;
	private int refPos;
	private	CigarOperator op;
	
	private char getRefBase()
		{
		if(this.genomicSequence==null)
			{
			return 'N';
			}
		else if(this.refPos>=1 && this.refPos<= this.genomicSequence.length())
			{
			return this.genomicSequence.charAt(this.refPos-1);
			}
		return '.';
		}
	
	private char getReadBase()
		{
		return this.readPos<0 || this.readbases==null  ||  this.readPos>=this.readbases.length?'.':(char)this.readbases[this.readPos];
		}
	private char getReadQual()
		{
		if(this.readPos<0 || this.readQuals==null || this.readQuals== SAMRecord.NULL_QUALS || this.readPos>=this.readQuals.length) return '.';
		return SAMUtils.phredToFastq(this.readQuals[this.readPos]);
		}
		
	private void writeAln()
		{
		this.out.print(this.rec.getReadName());
		this.out.print("\t");
		this.out.print(this.rec.getFlags());
		this.out.print("\t");
		this.out.print(this.rec.getMappingQuality());
		this.out.print("\t");
		this.out.print(this.rec.getReadUnmappedFlag()?SAMRecord.NO_ALIGNMENT_REFERENCE_NAME:this.rec.getReferenceName());
		this.out.print("\t");
		if(this.readPos!=-1)
			{
			final char c1 =this.getReadBase();
			this.out.print(this.readPos);
			this.out.print("\t");
			this.out.print(c1);
			this.out.print("\t");
			this.out.print(this.getReadQual());
			this.out.print("\t");
			}
		else
			{
			this.out.print(".\t.\t.\t");
			}
		
		if(this.refPos != -1)
			{
			char c3 = this.getRefBase();
			this.out.print(this.refPos);
			this.out.print("\t");
			this.out.print(c3);
			this.out.print("\t");
			}
		else
			{
			this.out.print(".\t.\t");
			}
		this.out.print(this.op==null?".":this.op.name());
		this.out.println();
		}
	
	private void printAln()
		{
		if(this.rec==null) return;
		
		
		this.readbases = rec.getReadBases();
		this.readQuals = rec.getBaseQualities();
		
		
		
		
		
		final Cigar cigar=rec.getCigar();
		if(cigar==null) return;


		if(this.indexedFastaSequenceFile!=null)
			{
			if(this.genomicSequence==null || !this.genomicSequence.getChrom().equals(rec.getContig()))
				{
				this.genomicSequence = new GenomicSequence(this.indexedFastaSequenceFile,rec.getContig());
				}
			}
		

		int readIndex = 0;
		int refIndex = rec.getUnclippedStart();
		 for (final CigarElement e : cigar.getCigarElements())
			 {
			 this.op = e.getOperator();
			 
			 switch (e.getOperator())
				 {
				 case S :
				 		{
				 		for(int i=0;i<e.getLength();++i)
				 			{
				 			this.readPos=  readIndex;
				 			this.refPos  = refIndex;
			 				writeAln();
				 			readIndex++;
				 			refIndex++;
				 			}
				 		break;
				 		}
				 case H :
					 	{
				 		for(int i=0;i<e.getLength();++i)
				 			{
				 			this.readPos=  -1;
				 			this.refPos  = refIndex;
			 				writeAln();
				 			refIndex++;
				 			}
						break; 
					 	}
				 case P : 
					 	{
					 	this.refPos  = -1;
					 	this.readPos = -1;
				 		for(int i=0;i<e.getLength();++i)
				 			{
				 			writeAln();
				 			}
						break; 
					 	}
				 case I :
				 		{
				 		this.refPos  = -1;
				 		for(int i=0;i<e.getLength();++i)
				 			{
				 			this.readPos=readIndex;
				 			writeAln();
				 			readIndex++;
				 			}
				 		break;
				 		}
				 case N :{
					 	if(this.skip_N) {
					 		refIndex += e.getLength();
					 		break;
					 		}
					 	else
					 		{
					 		this.readPos  = -1;
					 		for(int i=0;i<e.getLength();++i)
					 			{
					 			this.refPos = refIndex;
					 			writeAln();
					 			refIndex++;
					 			}
					 		}
					 	break;
				 		}
				 case D :
				 		{
				 		this.readPos  = -1;
				 		for(int i=0;i<e.getLength();++i)
				 			{
				 			this.refPos = refIndex;
				 			writeAln();
				 			refIndex++;
				 			}
				 		break;
				 		}
				 case M :
				 case EQ :
				 case X :
			 			{
				 		for(int i=0;i< e.getLength();++i)
				 			{
				 			this.refPos = refIndex;
				 			this.readPos = readIndex;
				 			writeAln();
				 			refIndex++;
				 			readIndex++;
				 			}
				 		break;
			 			}
				 default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
				 }
			 }
		}
	
	
	
	private void scan(final SAMFileHeader header,final CloseableIterator<SAMRecord> iter) 
		{
		long n_total=0L;
		long n_null_sequence  = 0L;
		long n_unmapped  = 0L;
		try{
			while(iter.hasNext())
				{
				this.rec = iter.next();
				
				if(this.rec.getReadUnmappedFlag()) {
					++n_unmapped;
					continue;
					}
				
				if(this.rec.getReadBases()==SAMRecord.NULL_SEQUENCE) {
					++n_null_sequence;
					continue;
					}
				
				printAln();
				++n_total;
				if(n_total%10_000==0 && this.out.checkError()) break;
				}
			if(n_null_sequence>0) LOG.warn("Ignored "+ n_null_sequence +" reads without sequence");
			if(n_unmapped>0) LOG.warn("Ignored "+ n_unmapped +" unmapped reads");
			}
		catch(final Exception err)
			{
			LOG.error("scan error:",err);
			throw new RuntimeException(String.valueOf(err.getMessage()),err);
			}
		finally
			{
			CloserUtil.close(iter);
			}
		}
	
	@Override
	protected int processInput(SAMFileHeader samHeader, CloseableIterator<SAMRecord> iter) {
		try
			{
			if(super.faidxPath!=null)
				{
				this.indexedFastaSequenceFile= ReferenceSequenceFileFactory.getReferenceSequenceFile(super.faidxPath);
				}
			
			
			this.out  =  openPathOrStdoutAsPrintWriter(this.outputFile);
			out.print("#Read-Name");
			out.print("\t");
			out.print("Flag");
			out.print("\t");
			out.print("MAPQ");
			out.print("\t");
			out.print("CHROM");
			out.print("\t");
			out.print("READ-POS0");
			out.print("\t");
			out.print("READ-BASE");
			out.print("\t");
			out.print("READ-QUAL");
			out.print("\t");
			out.print("REF-POS1");
			out.print("\t");
			out.print("REF-BASE");
			out.print("\t");
			out.print("CIGAR-OP");
			out.println();

			scan(samHeader,iter);			
			this.out.flush();this.out.close();this.out=null;
			return 0;
			}
		catch (final Throwable e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			CloserUtil.close(out);
			}
		}
	
	public static void main(final String[] args)
		{
		new Sam2Tsv().instanceMainWithExit(args);
		}
}
