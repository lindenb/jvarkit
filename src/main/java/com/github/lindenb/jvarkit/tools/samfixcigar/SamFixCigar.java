/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
* 2014 creation

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
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

@Program(name="samfixcigar",description="Fix Cigar String in SAM replacing 'M' by 'X' or '='")
public class SamFixCigar extends Launcher
	{
	private static final Logger LOG = Logger.build(SamFixCigar.class).make();

	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;


	@Parameter(names={"-r","--reference"},description="Indexed fasta Reference",required=true)
	private File faidx = null;
	
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();

	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	
	
	@Override
	public int doWork(List<String> args) {		
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
			LOG.info("Loading reference");
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(faidx);
			sfr = openSamReader(oneFileOrNull(args));
			final SAMFileHeader header=sfr.getFileHeader();
			sfw = this.writingBamArgs.openSAMFileWriter(outputFile,header, true);
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
					CigarOperator op = ce.getOperator();
					if(op.equals(CigarOperator.M))
						{
						for(int i=0;i< ce.getLength();++i)
    		    			{
							char c1=Character.toUpperCase((char)bases[readPos0]);
							char c2=Character.toUpperCase(refPos1-1< genomicSequence.length()?genomicSequence.charAt(refPos1-1):'*');
							
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
					CigarOperator op1 = newCigar.get(i).getOperator();
					int length1 = newCigar.get(i).getLength();
					
					if( i+1 <  newCigar.size() &&
						newCigar.get(i+1).getOperator()==op1)
						{
						CigarOperator op2= newCigar.get(i+1).getOperator();
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
		catch(Exception err)
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
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new SamFixCigar().instanceMainWithExit(args);

	}

}
