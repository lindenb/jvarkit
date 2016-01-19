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
package com.github.lindenb.jvarkit.tools.misc;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMValidationError;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

public class SamSlop extends AbstractSamSlop
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(SamSlop.class);
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		if(super.faidx==null)
			{
			return wrapException("Reference was not specified.");
			}
		if(super.defaultQual.length()!=1)
			{
			return wrapException("default quality should have length==1 "+super.defaultQual);
			}
		GenomicSequence genomicSequence=null;
		SamReader sfr=null;
		SAMFileWriter sfw=null;
		IndexedFastaSequenceFile indexedFastaSequenceFile=null;
		final char defaultQUAL=super.defaultQual.charAt(0);
		try
			{
			LOG.info("Loading reference");
			indexedFastaSequenceFile=new IndexedFastaSequenceFile(faidx);
			sfr = openSamReader(inputName);
			final SAMFileHeader header=sfr.getFileHeader();
			header.setSortOrder(SortOrder.unsorted);
			sfw = openSAMFileWriter(header, true);
			final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header);
			final SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				final SAMRecord rec=progress.watch(iter.next());
				final Cigar cigar=rec.getCigar();
				if( rec.getReadUnmappedFlag() ||
					cigar==null ||
					cigar.isEmpty() ||
					rec.getReadBases()==SAMRecord.NULL_SEQUENCE ||
					(super.extend5<=0 && super.extend3<=0)
					)
					{
					sfw.addAlignment(rec);
					continue;
					}
				
				final StringBuilder sbs = new StringBuilder(rec.getReadString());
				final StringBuilder sbq = new StringBuilder(rec.getBaseQualityString());

				
				if(genomicSequence==null ||
					genomicSequence.getSAMSequenceRecord().getSequenceIndex()!=rec.getReferenceIndex())
					{
					genomicSequence=new GenomicSequence(indexedFastaSequenceFile, rec.getReferenceName());
					}
				
				final int refPos1=(super.removeClip?rec.getAlignmentStart():rec.getUnclippedStart());
				final int endAlignmend1= (super.removeClip?rec.getAlignmentEnd():rec.getUnclippedEnd());
				final List<CigarElement> cl = new ArrayList<>(cigar.getCigarElements());
				
				if(!super.removeClip)
					{
					//replace clip S/H by match M
					for(int i=0;i< cl.size();++i)
						{
						final CigarElement ce = cl.get(i);
						if(ce.getOperator()==CigarOperator.S || ce.getOperator()==CigarOperator.H)
							{
							cl.set(i,new CigarElement(ce.getLength(), CigarOperator.M));
							}
						}
					}
				
				if(super.extend5>0 && refPos1>1)
					{
					
					if(super.removeClip)
						{
						///remove hard + soft clip 5'
						while(!cl.isEmpty())
							{
							//first
							final CigarElement ce = cl.get(0);
							//not a clip
							if(!(ce.getOperator()==CigarOperator.S || ce.getOperator()==CigarOperator.H))
								{
								break;
								}
							 if( ce.getOperator()==CigarOperator.S)
							 	{
								sbs.replace(0,ce.getLength(),"");
								sbq.replace(0,ce.getLength(),"");
							 	}
							cl.remove(0);//remove first
							}
						}
				
					final StringBuilder prefix= new StringBuilder(super.extend5);
					///append + soft clip 5'
					for(int i= 0;i< super.extend5;++i)
						{
						int x1 = (refPos1-1)-i;
						if(x1<1) break;//break if out of genome
						prefix.insert(0,genomicSequence.charAt(x1-1));
						}
					sbs.insert(0, prefix.toString());
					for(int i= 0;i<prefix.length();++i) sbq.insert(0, defaultQUAL);//preprend quality
					cl.add(0, new CigarElement(prefix.length(), CigarOperator.M));//prepend cigar
					rec.setAlignmentStart(refPos1-prefix.length());//update start pos
					}
				
				if(super.extend3>0 && rec.getAlignmentEnd()< genomicSequence.length())
					{
					if(super.removeClip)
						{
						///remove hard + soft clip 3'
						while(!cl.isEmpty())
							{
							//last
							final CigarElement ce = cl.get(cl.size()-1);
							//not a clip
							if(!(ce.getOperator()==CigarOperator.S || ce.getOperator()==CigarOperator.H))
								{
								break;
								}
							 if( ce.getOperator()==CigarOperator.S)
							 	{
								sbs.setLength(sbs.length()-ce.getLength());
								sbq.setLength(sbq.length()-ce.getLength());
							 	}
							 //remove last
							cl.remove(cl.size()-1);
							}
						}
					int extend=0;
					for(int pos1= endAlignmend1+1;
							pos1<= (endAlignmend1+super.extend3) && pos1<= genomicSequence.length() ;
							++pos1)
						{
						sbs.append(genomicSequence.charAt(pos1-1));
						sbq.append(defaultQUAL);
						++extend;
						}
					cl.add(new CigarElement(extend, CigarOperator.M));//append cigar
					}
				//simplify cigar
				int idx=0;
				while(idx+1<cl.size())
					{
					final CigarElement ce1 = cl.get(idx);
					final CigarElement ce2 = cl.get(idx+1);
					if(ce1.getOperator()==ce2.getOperator())
						{
						cl.set(idx,new CigarElement(ce1.getLength()+ce2.getLength(), ce1.getOperator()));
						cl.remove(idx+1);
						}
					else
						{
						idx++;
						}
					}
				
				
				rec.setCigar(new Cigar(cl));
				rec.setReadString(sbs.toString());
				rec.setBaseQualityString(sbq.toString());
				final List<SAMValidationError> errors = rec.isValid();
				if(errors!=null && !errors.isEmpty()) {
					for(SAMValidationError err:errors)
						{
						LOG.error(err.getMessage());
						}
				}
				
				//info("changed "+rec.getCigarString()+" to "+newCigarStr+" "+rec.getReadName()+" "+rec.getReadString());
				
				
				sfw.addAlignment(rec);
				}
			progress.finish();
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(indexedFastaSequenceFile);
			CloserUtil.close(sfr);
			CloserUtil.close(sfw);
			}
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new SamSlop().instanceMainWithExit(args);

	}

}
