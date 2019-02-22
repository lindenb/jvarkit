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


History:
* 2014-11 : handle clipped bases
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.sam2tsv;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
/**

BEGIN_DOC

### Output

Columns are:

 *  read name
 *  read flags
 *  reference name
 *  read-pos
 *  read-base
 *  read-qual
 *  ref-pos
 *  ref-base
 *  cigar-op



### Example
 


```
$ java -jar dist/sam2tsv.jar -A  \
    -r samtools-0.1.18/examples/toy.fa 
      samtools-0.1.18/examples/toy.sam
r001	163	ref	0	T	.	7	T	M
r001	163	ref	1	T	.	8	T	M
r001	163	ref	2	A	.	9	A	M
r001	163	ref	3	G	.	10	G	M
r001	163	ref	4	A	.	11	A	M
r001	163	ref	5	T	.	12	T	M
r001	163	ref	6	A	.	13	A	M
r001	163	ref	7	A	.	14	A	M
r001	163	ref	8	A	.	.	.	I
r001	163	ref	9	G	.	.	.	I
r001	163	ref	10	A	.	.	.	I
r001	163	ref	11	G	.	.	.	I
r001	163	ref	12	G	.	15	G	M
r001	163	ref	13	A	.	16	A	M
r001	163	ref	14	T	.	17	T	M
r001	163	ref	15	A	.	18	A	M
r001	163	ref	.	.	.	19	G	D
r001	163	ref	16	C	.	20	C	M
r001	163	ref	17	T	.	21	T	M
r001	163	ref	18	G	.	22	G	M
:   ref        7 TTAGATAAAGAGGATA-CTG 22      
:                ||||||||    |||| |||
:  r001        1 TTAGATAA----GATAGCTG 19      
r002	0	ref	1	A	.	.	.	I
r002	0	ref	2	A	.	.	.	I
r002	0	ref	3	A	.	9	A	M
r002	0	ref	4	G	.	10	G	M
r002	0	ref	5	A	.	11	A	M
r002	0	ref	6	T	.	12	T	M
r002	0	ref	7	A	.	13	A	M
r002	0	ref	8	A	.	14	A	M
r002	0	ref	9	G	.	.	.	I
r002	0	ref	10	G	.	.	.	I
r002	0	ref	11	G	.	15	G	M
r002	0	ref	12	A	.	16	A	M
r002	0	ref	13	T	.	17	T	M  
(...)   

```


## Example 2

sam2tsv can read data from a linux pipe.

```
samtools view -h input.bam | java -jar dist/sam2tsv.jar
```




### History

 *  Moved to a standard argc/argv command line
 *  2014-04: added qual and samflag. Fixed a bug in soft-clip
 *  2014-11: manage hard+soft clip
 *  2019-02 : manage reads without qualities, contig name converter

### Citations


Sam2tsv was cited in : 

  * "Illumina TruSeq Synthetic Long-Reads Empower De Novo Assembly and Resolve Complex, Highly-Repetitive Transposable Elements" . McCoy RC, Taylor RW, Blauwkamp TA, Kelley JL, Kertesz M, et al. (2014) Illumina TruSeq Synthetic Long-Reads Empower De Novo Assembly and Resolve Complex, Highly-Repetitive Transposable Elements. PLoS ONE 9(9): e106689. doi: 10.1371/journal.pone.0106689  http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0106689
  * "High-Throughput Identification of Genetic Variation Impact on pre-mRNA Splicing Efficiency". Scott I Adamson, Lijun Zhan, Brenton R Graveley. doi: [https://doi.org/10.1101/191122](https://doi.org/10.1101/191122).
  * "Linkage of A-to-I RNA editing in metazoans and the impact on genome evolution "  Molecular Biology and Evolution, msx274, https://doi.org/10.1093/molbev/msx274
  * "Vex-seq: high-throughput identification of the impact of genetic variation on pre-mRNA splicing efficiency" Genome Biology201819:71 https://doi.org/10.1186/s13059-018-1437-x
  * "Accurate detection of m6A RNA modifications in native RNA sequences" Huanle Liu, Oguzhan Begik, Morghan C Lucas, Christopher E Mason, Schraga Schwartz, John S Mattick, Martin A Smith, Eva Maria Novoa bioRxiv 525741; doi: https://doi.org/10.1101/525741 

END_DOC
*/
@Program(name="sam2tsv",
	description="Prints the SAM alignments as a TAB delimited file.",
	keywords={"sam","bam","table","tsv"},
	biostars={157232,59647,253828,264875,277493},
	modificationDate="20190222")
public class Sam2Tsv
	extends Launcher
	{
	private static final Logger LOG = Logger.build(Sam2Tsv.class).make();

	@Parameter(names={"-o","--output"},description= OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-A","--printAlignments"},description="Print Alignments")
	private boolean printAlignment = false;
	@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private File refFile = null;

	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private GenomicSequence genomicSequence=null;
	private SAMSequenceDictionary refDict = null;
	private ContigNameConverter contigNameConverter = null;
	/** lines for alignments */
	private StringBuilder L1=null;
	private StringBuilder L2=null;
	private StringBuilder L3=null;

	private PrintWriter out = null;
	
	private class Row
		{
		SAMRecord rec;
		byte readbases[];
		byte readQuals[];

		int readPos;
		int refPos;
		CigarOperator op;
		
		public char getRefBase()
			{
			if(Sam2Tsv.this.genomicSequence==null)
				{
				return 'N';
				}
			else if(refPos>=1 && refPos<= genomicSequence.length())
 				{
				return genomicSequence.charAt(refPos-1);
 				}
			return '.';
			}
		
		public char getReadBase()
			{
			return readPos==-1 || this.readPos>=this.readbases.length?'.':(char)this.readbases[this.readPos];
			}
		public char getReadQual()
			{
			byte c= readPos==-1 || this.readQuals==null || this.readPos>=this.readQuals.length?(byte)0:this.readQuals[this.readPos];
			return SAMUtils.phredToFastq(c);
			}
		}
	
	private void writeAln(final Row row)
			{
			char c1;
			char c3;
			this.out.print(row.rec.getReadName());
			this.out.print("\t");
			this.out.print(row.rec.getFlags());
			this.out.print("\t");
			this.out.print(row.rec.getReadUnmappedFlag()?".":row.rec.getReferenceName());
			this.out.print("\t");
			if(row.readPos!=-1)
				{
				c1 =row.getReadBase();
				this.out.print(row.readPos);
				this.out.print("\t");
				this.out.print(c1);
				this.out.print("\t");
				this.out.print(row.getReadQual());
				this.out.print("\t");
				}
			else
				{
				c1= '-';
				this.out.print(".\t.\t.\t");
				}
			
			if(row.refPos != -1)
				{
				c3 = row.getRefBase();
				this.out.print(row.refPos);
				this.out.print("\t");
				this.out.print(c3);
				this.out.print("\t");
				}
			else
				{
				c3= '-';
				this.out.print(".\t.\t");
				}
			this.out.print(row.op==null?".":row.op.name());
			this.out.println();
			
			if(this.printAlignment)
				{
				L1.append(c1);
				L3.append(c3);
				
				if(Character.isLetter(c1) &&  Character.toUpperCase(c1)== Character.toUpperCase(c3))
					{
					L2.append('|');
					}
				else
					{
					L2.append(' ');
					}
				}
			}
	
	private void printAln(final Row row)
		{
		final SAMRecord rec = row.rec;
		if(rec==null) return;
		final Cigar cigar=rec.getCigar();
		if(cigar==null) return;
		
		row.readbases = rec.getReadBases();
		row.readQuals = rec.getBaseQualities()==SAMRecord.NULL_QUALS?
				StringUtils.repeat(row.readbases.length,'#').getBytes():
				rec.getBaseQualities()
				;
		if(row.readbases==null )
			{
			row.op=null;
			row.refPos=-1;
			row.readPos=-1;
			writeAln(row);
			return;
			}
		if(rec.getReadUnmappedFlag())
			{
			row.op=null;
			row.refPos=-1;
			for(int i=0;i< row.readbases.length;++i)
				{
				row.readPos=i;
				writeAln(row);
				}
			return;
			}
		
		//fix hard clipped reads
		final StringBuilder fixReadBases=new StringBuilder(row.readbases.length);
		final  StringBuilder fixReadQuals=new StringBuilder(row.readbases.length);
		int readIndex = 0;
		for (final CigarElement ce : cigar.getCigarElements())
			 {
			 final CigarOperator op= ce.getOperator();
			 
			 for(int i=0;i< ce.getLength();++i)
				{
				if(op.equals(CigarOperator.H))
					{
					fixReadBases.append('*');
					fixReadQuals.append('*');
					}
				else if(!op.consumesReadBases())
					{
					break;
					}
				else
					{
					fixReadBases.append((char)row.readbases[readIndex]);
					fixReadQuals.append(
							row.readQuals==null ||
							row.readQuals.length<=readIndex ?
							'*':(char)row.readQuals[readIndex]);
					readIndex++;
					}
				}
			 }
		row.readbases = fixReadBases.toString().getBytes();
		row.readQuals = fixReadQuals.toString().getBytes();

		if(this.indexedFastaSequenceFile!=null)
			{
			final String ctg = this.contigNameConverter.apply(rec.getContig());
			if(StringUtils.isBlank(ctg)) throw new JvarkitException.ContigNotFoundInDictionary(rec.getContig(),this.refDict);
			if(this.genomicSequence==null || !this.genomicSequence.getChrom().equals(ctg))
				{
				this.genomicSequence = new GenomicSequence(this.indexedFastaSequenceFile,ctg);
				}
			}
		

		 readIndex = 0;
		 int refIndex = rec.getUnclippedStart();
		 				 
		 for (final CigarElement e : cigar.getCigarElements())
			 {
			 row.op=e.getOperator();
			 
			 switch (e.getOperator())
				 {
				 case S :
				 case H : //length of read has been fixed previously, so same as 'S'
					 	{
					 	
				 		for(int i=0;i<e.getLength();++i)
				 			{
				 			row.readPos=  readIndex;
				 			row.refPos  = refIndex;
			 				writeAln(row);
				 			readIndex++;
				 			refIndex++;//because we used getUnclippedStart
				 			}
						break; 
					 	}
				 case P : 
					 	{
					 	row.refPos  = -1;
					 	row.readPos = -1;
				 		for(int i=0;i<e.getLength();++i)
				 			{
				 			writeAln(row);
				 			}
						break; 
					 	}
				 case I :
				 		{
				 		row.refPos  = -1;
				 		for(int i=0;i<e.getLength();++i)
				 			{
				 			row.readPos=readIndex;
				 			writeAln(row);
				 			readIndex++;
				 			}
				 		break;
				 		}
				 case N :  //cont. -- reference skip
				 case D :
				 		{
				 		row.readPos  = -1;
				 		for(int i=0;i<e.getLength();++i)
				 			{
				 			row.refPos = refIndex;
				 			writeAln(row);
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
				 			row.refPos = refIndex;
				 			row.readPos = readIndex;
				 			writeAln(row);
				 			refIndex++;
				 			readIndex++;
				 			}
				 		break;
			 			}
					
				 default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
				 }

			 }
	
		
		
		 if(printAlignment)
				{
				
				final int len=Math.max(rec.getReadNameLength(), rec.getReferenceName().length())+2;

				this.out.printf(":%"+len+"s %8d %s %-8d\n",
						rec.getReferenceName(),
						rec.getUnclippedStart(),
						L3.toString(),
						rec.getUnclippedEnd()
						);
				this.out.printf(":%"+len+"s %8s %s\n",
						"",
						"",
						L2.toString()
						);

				this.out.printf(":%"+len+"s %8d %s %-8d\n",
						rec.getReadName(),
						1,
						L1.toString(),
						rec.getReadLength()
						);

				L1.setLength(0);
				L2.setLength(0);
				L3.setLength(0);
				}

		}
	
	
	
	private void scan(final SamReader r) 
		{
		final Row row=new Row();
		SAMRecordIterator iter=null;
		try{
			final ProgressFactory.Watcher<SAMRecord> progress= ProgressFactory.newInstance().dictionary(r.getFileHeader()).logger(LOG).build();
			iter=r.iterator();	
			while(iter.hasNext())
				{
				row.rec =progress.apply(iter.next());
				if(row.rec.getReadBases()==SAMRecord.NULL_SEQUENCE) {
					LOG.warn("Ignoring read without sequence: "+row.rec.getReadName());
					continue;
					}
				printAln(row);
				if(this.out.checkError()) break;
				}
			progress.close();
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
	public int doWork(final List<String> args) {
		
		if(printAlignment)
			{
			L1=new StringBuilder();
			L2=new StringBuilder();
			L3=new StringBuilder();
			}
		
		SamReader samFileReader=null;
		try
			{
			if(this.refFile!=null)
				{
				this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(refFile);
				this.refDict = SequenceDictionaryUtils.extractRequired(this.indexedFastaSequenceFile);
				this.contigNameConverter = ContigNameConverter.fromOneDictionary(this.refDict);
				}
			this.out  =  openFileOrStdoutAsPrintWriter(outputFile);
			this.out.println("#READ_NAME\tFLAG\tCHROM\tREAD_POS\tBASE\tQUAL\tREF_POS\tREF\tOP");
			samFileReader= openSamReader(oneFileOrNull(args));
			
			scan(samFileReader);
			samFileReader.close();
			samFileReader = null;
			this.out.flush();this.out.close();this.out=null;
			return RETURN_OK;
			}
		catch (final Throwable e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			CloserUtil.close(samFileReader);
			CloserUtil.close(out);
			L1=null;
			L2=null;
			L3=null;
			}
		}
	
	public static void main(final String[] args)
		{
		new Sam2Tsv().instanceMainWithExit(args);
		}
}
