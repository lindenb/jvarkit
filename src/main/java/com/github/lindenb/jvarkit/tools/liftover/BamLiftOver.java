/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.liftover;

import java.util.Collection;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;

public class BamLiftOver extends AbstractBamLiftOver
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(BamLiftOver.class);

	@Override
	protected Collection<Throwable> call(final String inputName) throws Exception {
		final double minMatch=(super.userMinMatch<=0.0?LiftOver.DEFAULT_LIFTOVER_MINMATCH:super.userMinMatch);
		if(super.liftOverFile==null)
			{
			return wrapException("LiftOver file is undefined.");
			}
		if(super.faidx==null)
			{
			return wrapException("New Sequence Dictionary file is undefined.");
			}
		SAMRecordIterator iter=null;
		SamReader sfr=null;
		SAMFileWriter sfw=null;
		try
			{
			LOG.info("Reading "+liftOverFile);
			LiftOver liftOver=new LiftOver(liftOverFile);
			liftOver.setLiftOverMinMatch(minMatch);

			
			final SAMSequenceDictionary newDict=new SAMSequenceDictionaryFactory().load(super.faidx);
			
			
			sfr=super.openSamReader(inputName);
			

			final SAMFileHeader headerIn=sfr.getFileHeader();
			final SAMFileHeader headerOut=headerIn.clone();
			headerOut.setSortOrder(SortOrder.unsorted);
			
			sfw = openSAMFileWriter(headerOut, true);
			
			
			iter=sfr.iterator();
			while(iter.hasNext())
				{
				final SAMRecord rec=iter.next();
				final SAMRecord copy=(SAMRecord)rec.clone();
				copy.setHeader(headerOut);
				final StringBuilder sb=new StringBuilder();
				if(!rec.getReadUnmappedFlag())
					{
					final String chrom=rec.getReferenceName();
					int pos=rec.getAlignmentStart();
					final Interval interval=liftOver.liftOver(new Interval(chrom, pos,pos,rec.getReadNegativeStrandFlag(),null));
					if(interval!=null)
						{
						sb.append(chrom+":"+pos+":"+(rec.getReadNegativeStrandFlag()?"-":"+"));
						final SAMSequenceRecord ssr=newDict.getSequence(interval.getContig());
						if(ssr==null)
							{
							sfr.close();
							sfr=null;
							return wrapException("the chromosome "+interval.getContig()+" is undefined in the sequence dict.");
							}
						copy.setReferenceName(ssr.getSequenceName());
						copy.setReferenceIndex(ssr.getSequenceIndex());
						copy.setAlignmentStart(interval.getStart());
						copy.setReadNegativeStrandFlag(interval.isNegativeStrand());
						if(rec.getReadNegativeStrandFlag()!=copy.getReadNegativeStrandFlag()) {
							copy.setReadString(AcidNucleics.reverseComplement(rec.getReadString()));
							
							byte qual[]= rec.getBaseQualities();
							byte quals2[]=  new byte[qual.length];
							for(int i=0;i< qual.length;++i) {
								quals2[i]=qual[(qual.length-1)-i];
							}
							copy.setBaseQualities(quals2);
							}
						}
					else
						{
						sb.append(".");
						SAMUtils.makeReadUnmapped(copy);
						}
					}
				
				
				if(rec.getReadPairedFlag() && !rec.getMateUnmappedFlag())
					{
					sb.append("/");
					String chrom=rec.getMateReferenceName();
					int pos=rec.getMateAlignmentStart();
					final Interval interval=liftOver.liftOver(new Interval(chrom, pos,pos,rec.getMateNegativeStrandFlag(),null));
					if(interval!=null)
						{
						sb.append(chrom+":"+pos+":"+(rec.getMateNegativeStrandFlag()?"-":"+"));
						final SAMSequenceRecord ssr=newDict.getSequence(interval.getContig());
						if(ssr==null)
							{
							sfr.close();
							sfr=null;
							return wrapException("the chromosome "+interval.getContig()+" is undefined in the sequence dict.");
							}
						copy.setMateReferenceName(ssr.getSequenceName());
						copy.setMateReferenceIndex(ssr.getSequenceIndex());
						copy.setMateAlignmentStart(interval.getStart());
						copy.setMateNegativeStrandFlag(interval.isNegativeStrand());
						
						if(!copy.getReadUnmappedFlag() &&
							copy.getReferenceIndex()==copy.getMateReferenceIndex() 
							// && copy.getReadNegativeStrandFlag()!=copy.getMateNegativeStrandFlag()
							)
							{
							//don't change ?
							}
						else
							{
							copy.setProperPairFlag(false);
							copy.setInferredInsertSize(0);
							}
						}
					else
						{
						sb.append(".");
						SAMUtils.makeReadUnmapped(copy);
						}
					}
				if(sb.length()>0) copy.setAttribute("LO", sb.toString());
				sfw.addAlignment(copy);
				}
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(sfr);
			CloserUtil.close(sfw);
			}
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new BamLiftOver().instanceMainWithExit(args);
		}

	}
