/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.sam2tsv;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.function.UnaryOperator;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.MultiBamLauncher;
import com.github.lindenb.jvarkit.samtools.SAMRecordLeftAligner;
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
  * Yang & al. "Sequencing 5-Formyluracil in Genomic DNA at Single-Base Resolution"  (2021) Analytical Chemistry doi: 10.1021/acs.analchem.1c03339
  * Xiao Shu, Chenyang Huang, Tengwei Li, Jie Cao, Jianzhao Liu,.a6A-seq: N6-allyladenosine-based cellular messenger RNA metabolic labelling and sequencing,.Fundamental Research,.2023,.,.ISSN 2667-3258, https://doi.org/10.1016/j.fmre.2023.04.010..(https://www.sciencedirect.com/science/article/pii/S2667325823001267).
  * Xu, Z., Sziraki, A., Lee, J. et al. Dissecting key regulators of transcriptome kinetics through scalable single-cell RNA profiling of pooled CRISPR screens. Nat Biotechnol (2023). https://doi.org/10.1038/s41587-023-01948-9
  * Braton & al (2023) HIVepsilon-seq-scalable characterization of intact persistent proviral HIV reservoirs in women. Virology. https://doi.org/10.1128/jvi.00705-23
  * Yin, K., Xu, Y., Guo, Y. et al. Dyna-vivo-seq unveils cellular RNA dynamics during acute kidney injury via in vivo metabolic RNA labeling-based scRNA-seq. Nat Commun 15, 9866 (2024). https://doi.org/10.1038/s41467-024-54202-4
  * Heng, J. J. W. (2024). Development of sequencing methods to illuminate the epitranscriptome. Doctoral thesis, Nanyang Technological University, Singapore. https://hdl.handle.net/10356/182139
  * Liu, C., Liang, H., Wan, A.H. et al. Decoding the m6A epitranscriptomic landscape for biotechnological applications using a direct RNA sequencing approach. Nat Commun 16, 798 (2025). https://doi.org/10.1038/s41467-025-56173-6

END_DOC
*/
@Program(name="sam2tsv",
	description="Prints the SAM alignments as a TAB delimited file.",
	keywords={"sam","bam","table","cram","tsv"},
	biostars={157232,59647,253828,264875,277493},
	modificationDate="20250327",
	creationDate="20170712",
	jvarkit_amalgamion =  true,
	nfcore="https://nf-co.re/modules/jvarkit_sam2tsv.html",
	menu="BAM Manipulation"
	)
public class Sam2Tsv
	extends MultiBamLauncher
	{
	private static final Logger LOG = Logger.of(Sam2Tsv.class);
	private static final int NO_POSITION=-1;
	@Parameter(names={"-o","--output"},description= OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-N","--skip-N"},description="Skip 'N' operator")
	private boolean skip_N=false;
	@Parameter(names={"--left-align"},description=SAMRecordLeftAligner.OPT_DESC)
	private boolean enable_left_alin=false;


	private ReferenceSequenceFile indexedFastaSequenceFile=null;
	private GenomicSequence genomicSequence=null;
	private UnaryOperator<SAMRecord> leftAligner= R->R;
	
	private char getRefBase(final int refPos)
		{
		if(this.genomicSequence==null)
			{
			return 'N';
			}
		else if(refPos>=1 &&  refPos<= this.genomicSequence.length())
			{
			return this.genomicSequence.charAt(refPos-1);
			}
		return '.';
		}
	
	private char getReadBase(final SAMRecord rec,int readPos,byte[] readbases)
		{
		if(readbases==null  || readbases==SAMRecord.NULL_SEQUENCE ) return '.';
		if( readPos<0 || readPos>=readbases.length) {
			System.err.println("Error with getReadBase "+rec.getPairedReadName()+" readpos="+readPos+" and readbases.length="+readbases.length);
			return '.';
			}
		return (char)readbases[readPos];
		}
	private char getReadQual(int readPos,byte[] readQuals)
		{
		if(readPos<0 || readQuals==null ||readQuals== SAMRecord.NULL_QUALS || readPos>=readQuals.length) return '.';
		return SAMUtils.phredToFastq(readQuals[readPos]);
		}
		
	private void writeAln(
			final PrintWriter pw,
			final SAMRecord rec, 
			final CigarOperator op,
			final byte[] readbases,
			final byte[] readQuals,
			final int readPos,
			final int refPos
			)
		{
		pw.print(rec.getReadName());
		pw.print("\t");
		pw.print(rec.getFlags());
		pw.print("\t");
		pw.print(rec.getMappingQuality());
		pw.print("\t");
		pw.print(rec.getReadUnmappedFlag()?SAMRecord.NO_ALIGNMENT_REFERENCE_NAME:rec.getReferenceName());
		pw.print("\t");
		if(readPos!= NO_POSITION )
			{
			final char c1 =this.getReadBase(rec,readPos,readbases);
			pw.print(readPos);
			pw.print("\t");
			pw.print(c1);
			pw.print("\t");
			pw.print(this.getReadQual(readPos,readQuals));
			pw.print("\t");
			}
		else
			{
			pw.print(".\t.\t.\t");
			}
		
		if(refPos != NO_POSITION)
			{
			char c3 = this.getRefBase(refPos);
			pw.print(refPos);
			pw.print("\t");
			pw.print(c3);
			pw.print("\t");
			}
		else
			{
			pw.print(".\t.\t");
			}
		pw.print(op==null?".":op.name());
		pw.println();
		}
	
	private void printAln(final PrintWriter pw,final SAMRecord rec)
		{
		if(rec==null) return;
		
		
		final byte[] readbases = rec.getReadBases();
		final byte[] readQuals = rec.getBaseQualities();
		
		
		final Cigar cigar=rec.getCigar();
		if(cigar==null) return;


		if(this.indexedFastaSequenceFile!=null)
			{
			if(this.genomicSequence==null || !this.genomicSequence.hasName(rec.getContig()))
				{
				this.genomicSequence = new GenomicSequence(this.indexedFastaSequenceFile,rec.getContig());
				}
			}
		

		int readIndex = 0;
		int refIndex = rec.getUnclippedStart();
		for (final CigarElement e : cigar.getCigarElements())
		 {
		final  CigarOperator op = e.getOperator();
		 
		 switch (op)
			 {
			 case H :
				 	{
			 		for(int i=0;i<e.getLength();++i)
			 			{
		 				writeAln(pw,rec,op,readbases,readQuals,NO_POSITION,refIndex);
			 			refIndex++;
			 			}
					break; 
				 	}
			 case P : 
				 	{
			 		for(int i=0;i<e.getLength();++i)
			 			{
			 			writeAln(pw,rec,op,readbases,readQuals,NO_POSITION,NO_POSITION);
			 			}
					break; 
				 	}
			 case I :
			 		{
			 		for(int i=0;i<e.getLength();++i)
			 			{
			 			writeAln(pw,rec,op,readbases,readQuals,readIndex,NO_POSITION);
			 			readIndex++;
			 			}
			 		break;
			 		}
			 case N :
			 case D :
			 		{
		 			if(this.skip_N && op.equals(CigarOperator.N)) {
				 		refIndex += e.getLength();
				 		}
		 			else
			 			{
				 		for(int i=0;i<e.getLength();++i)
				 			{
				 			writeAln(pw,rec,op,readbases,readQuals,NO_POSITION,refIndex);
				 			refIndex++;
				 			}
			 			}
			 		break;
			 		}
			 case S:
			 case M :
			 case EQ :
			 case X :
		 			{
			 		for(int i=0;i< e.getLength();++i)
			 			{
			 			writeAln(pw,rec,op,readbases,readQuals,readIndex,refIndex);
			 			refIndex++;
			 			readIndex++;
			 			}
			 		break;
		 			}
			 default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
			 }
		 }
		}
	
	@Override
	protected int beforeSam() {
		try {
			if(super.faidxPath!=null)
				{
				this.indexedFastaSequenceFile= ReferenceSequenceFileFactory.getReferenceSequenceFile(super.faidxPath);
				}
			}
		catch(Throwable err) {
			getLogger().error(err);
			return -1;
			}
		if(this.enable_left_alin) {
			if(this.indexedFastaSequenceFile==null) {
				getLogger().error("Reference is required with using --left-align");
				return -1;
				}
			this.leftAligner = new SAMRecordLeftAligner(indexedFastaSequenceFile);
			}
		else
			{
			this.leftAligner =  R->R;
			}
		return super.beforeSam();
		}
	
	
	@Override
	protected void afterSam() {
		if(leftAligner!=null) {
			leftAligner=null;
			}
		if(this.indexedFastaSequenceFile!=null) {
			try { this.indexedFastaSequenceFile.close();}
			catch(Throwable err) {}
			this.indexedFastaSequenceFile=null;
			}
		super.afterSam();
		}
	
	@Override
	protected int processInput(final SAMFileHeader samHeader, final CloseableIterator<SAMRecord> iter) {
		try
			{
			try(PrintWriter out  =  openPathOrStdoutAsPrintWriter(this.outputFile)) {
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
				
				
				
				long n_total=0L;
				long n_null_sequence  = 0L;
				long n_unmapped  = 0L;
				while(iter.hasNext())
					{
					final SAMRecord rec = this.leftAligner.apply(iter.next());
					
					if(rec.getReadUnmappedFlag()) {
						++n_unmapped;
						continue;
						}
					
					if(rec.getReadBases()==SAMRecord.NULL_SEQUENCE) {
						++n_null_sequence;
						continue;
						}
					
					printAln(out,rec);
					++n_total;
					if(n_total%10_000==0 && out.checkError()) break;
					}
				if(n_null_sequence>0) LOG.warn("Ignored "+ n_null_sequence +" reads without sequence");
				if(n_unmapped>0) LOG.warn("Ignored "+ n_unmapped +" unmapped reads");
				
				out.flush();
				}
			return 0;
			}
		catch (final Throwable e)
			{
			LOG.error(e);
			return -1;
			}
		}
	
	public static void main(final String[] args)
		{
		new Sam2Tsv().instanceMainWithExit(args);
		}
}
