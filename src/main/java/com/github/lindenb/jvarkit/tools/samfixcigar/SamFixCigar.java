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
package com.github.lindenb.jvarkit.tools.samfixcigar;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
/**

BEGIN_DOC


### Example


the input file


```
$ cat toy.sam

@SQ     SN:ref  LN:45
@SQ     SN:ref2 LN:40
r001    163     ref     7       30      8M4I4M1D3M      =       37      39      TTAGATAAAGAGGATACTG     *       XX:B:S,12561,2,20,112
r002    0       ref     9       30      1S2I6M1P1I1P1I4M2I      *       0       0       AAAAGATAAGGGATAAA       *
r003    0       ref     9       30      5H6M    *       0       0       AGCTAA  *
r004    0       ref     16      30      6M14N1I5M       *       0       0       ATAGCTCTCAGC    *
r003    16      ref     29      30      6H5M    *       0       0       TAGGC   *
r001    83      ref     37      30      9M      =       7       -39     CAGCGCCAT       *
x1      0       ref2    1       30      20M     *       0       0       aggttttataaaacaaataa    ????????????????????
x2      0       ref2    2       30      21M     *       0       0       ggttttataaaacaaataatt   ?????????????????????
x3      0       ref2    6       30      9M4I13M *       0       0       ttataaaacAAATaattaagtctaca      ??????????????????????????
x4      0       ref2    10      30      25M     *       0       0       CaaaTaattaagtctacagagcaac       ?????????????????????????
x5      0       ref2    12      30      24M     *       0       0       aaTaattaagtctacagagcaact        ????????????????????????
x6      0       ref2    14      30      23M     *       0       0       Taattaagtctacagagcaacta ???????????????????????

```


processing with samfixcigar


```
$ java -jar dist/samfixcigar.jar \
     -r samtools-0.1.19/examples/toy.fa \
     samtools-0.1.19/examples/toy.sam
@HD     VN:1.4  SO:unsorted
@SQ     SN:ref  LN:45
@SQ     SN:ref2 LN:40
r001    163     ref     7       30      8=4I4=1D3=      =       37      39      TTAGATAAAGAGGATACTG     *       XX:B:S,12561,2,20,112
r002    0       ref     9       30      1S2I6=1P1I1P1I1X1=2X2I  *       0       0       AAAAGATAAGGGATAAA       *
r003    0       ref     9       30      2=1X3=  *       0       0       AGCTAA  *
r004    0       ref     16      30      6=14N1I5=       *       0       0       ATAGCTCTCAGC    *
r003    16      ref     29      30      5=      *       0       0       TAGGC   *
r001    83      ref     37      30      9=      =       7       -39     CAGCGCCAT       *
x1      0       ref2    1       30      16=1X3= *       0       0       AGGTTTTATAAAACAAATAA    ????????????????????
x2      0       ref2    2       30      15=1X3=1X1=     *       0       0       GGTTTTATAAAACAAATAATT   ?????????????????????
x3      0       ref2    6       30      9=4I13= *       0       0       TTATAAAACAAATAATTAAGTCTACA      ??????????????????????????
x4      0       ref2    10      30      1X3=1X20=       *       0       0       CAAATAATTAAGTCTACAGAGCAAC       ?????????????????????????
x5      0       ref2    12      30      2=1X21= *       0       0       AATAATTAAGTCTACAGAGCAACT        ????????????????????????
x6      0       ref2    14      30      1X22=   *       0       0       TAATTAAGTCTACAGAGCAACTA ???????????????????????
```

### Usage in the literature

This tool was cited in Extensive sequencing of seven human genomes to characterize benchmark reference materials Sci Data. 2016; 3: 160025..

END_DOC

*/

@Program(name="samfixcigar",
	description="Fix Cigar String in SAM replacing 'M' by 'X' or '='",
	keywords={"sam","bam","cigar"},
	biostars= {312430,340479}
	)
public class SamFixCigar extends Launcher
	{
	private static final Logger LOG = Logger.build(SamFixCigar.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private File faidx = null;
	
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();

	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	
	
	@Override
	public int doWork(final List<String> args) {		
		if(this.faidx==null)
			{
			LOG.error("Reference was not specified.");
			return -1;
			}
		GenomicSequence genomicSequence=null;

		SamReader sfr=null;
		SAMFileWriter sfw=null;
		try
			{
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(faidx);
			sfr = openSamReader(oneFileOrNull(args));
			final SAMFileHeader header = sfr.getFileHeader();
			header.addComment("Processed with "+getProgramName()+" "+JVarkitVersion.getInstance().getLabel());
			sfw = this.writingBamArgs.
					setReferenceFile(this.faidx).
					openSAMFileWriter(outputFile,header, true);
			final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header);
			final List<CigarElement> newCigar=new ArrayList<CigarElement>();
			final SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				final SAMRecord rec=progress.watch(iter.next());
				Cigar cigar=rec.getCigar();
				byte bases[]=rec.getReadBases();
				if( rec.getReadUnmappedFlag() ||
					cigar==null ||
					cigar.getCigarElements().isEmpty() ||
					bases==null)
					{
					sfw.addAlignment(rec);
					continue;
					}
				
				if(genomicSequence==null ||
					genomicSequence.getSAMSequenceRecord().getSequenceIndex()!=rec.getReferenceIndex())
					{
					genomicSequence=new GenomicSequence(indexedFastaSequenceFile, rec.getReferenceName());
					}
				
				newCigar.clear();
				int refPos1=rec.getAlignmentStart();
				int readPos0=0;
				
				for(final CigarElement ce:cigar.getCigarElements())
					{
					final CigarOperator op = ce.getOperator();
					if(op.equals(CigarOperator.M))
						{
						for(int i=0;i< ce.getLength();++i)
    		    			{
							final char c1=Character.toUpperCase((char)bases[readPos0]);
							final char c2=Character.toUpperCase(refPos1-1< genomicSequence.length()?genomicSequence.charAt(refPos1-1):'*');
							
							if(c2=='N' || c1==c2)
								{
								newCigar.add(new CigarElement(1, CigarOperator.EQ));
								}
							else
								{
								newCigar.add(new CigarElement(1, CigarOperator.X));
								}
    						refPos1++;
    						readPos0++;
		    				}
						}
					else
						{
						newCigar.add(ce);
						if(op.consumesReadBases()) readPos0+=ce.getLength();	
						if(op.consumesReferenceBases()) refPos1+=ce.getLength();	
						}
					}
				
				int i=0;
				while(i< newCigar.size())
					{
					final CigarOperator op1 = newCigar.get(i).getOperator();
					final int length1 = newCigar.get(i).getLength();
					
					if( i+1 <  newCigar.size() &&
						newCigar.get(i+1).getOperator()==op1)
						{
						final CigarOperator op2= newCigar.get(i+1).getOperator();
						int length2=newCigar.get(i+1).getLength();

						 newCigar.set(i,new CigarElement(length1+length2, op2));
						 newCigar.remove(i+1);
						}
					else
						{
						++i;
						}
					}
				cigar=new Cigar(newCigar);
				//info("changed "+rec.getCigarString()+" to "+newCigarStr+" "+rec.getReadName()+" "+rec.getReadString());
				rec.setCigar(cigar);
				
				sfw.addAlignment(rec);
				}
			progress.finish();
			return RETURN_OK;
			}
		catch(final  Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			CloserUtil.close(sfr);
			CloserUtil.close(sfw);
			}
		}
	
	public static void main(final  String[] args) {
		new SamFixCigar().instanceMainWithExit(args);

	}

}
