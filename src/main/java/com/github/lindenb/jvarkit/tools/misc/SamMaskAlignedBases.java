
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


History:
* 2016 : creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.util.Collection;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.CloserUtil;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class SamMaskAlignedBases
	extends AbstractSamMaskAlignedBases
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(SamMaskAlignedBases.class);
	
	
	public SamMaskAlignedBases()
	{
	}
	

	@Override
	public Collection<Throwable> call(final String inputName) throws Exception 
		{
		final byte RESET_CHAR=(byte)'N';
		final byte RESET_QUAL=(byte)SAMUtils.fastqToPhred('#');
		
		long nRecords=0L;
		long nRecordMasked=0L;
		long nBasesMasked=0L;
		long nBases=0L;
		SAMRecordIterator iter=null;
		SamReader sfr=null;
		SAMFileWriter sfw =null;
		try
			{			
			sfr = openSamReader(inputName);
			
			final SAMFileHeader header1=sfr.getFileHeader();
			
			if(header1==null)
				{
				return wrapException("File header missing");
				}
			final SAMFileHeader header2=header1.clone();
			header2.addComment(getName()+":"+getVersion()+":"+getProgramCommandLine());
			header2.setSortOrder(SortOrder.unsorted);
			
			sfw =  openSAMFileWriter(header2, true);
			
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header1);
			iter = sfr.iterator();
			while(iter.hasNext())
				{
				final SAMRecord rec=progress.watch(iter.next());
				++nRecords;
				nBases+=rec.getReadLength();
				if(rec.getReadUnmappedFlag()) {
					SAMUtils.makeReadUnmapped(rec);
					sfw.addAlignment(rec);
					continue;
				}
				
				if(rec.isSecondaryOrSupplementary()) {
					continue;
				}
				
				
				final Cigar cigar = rec.getCigar();
				byte bases[] = rec.getReadBases();
				byte quals[] = rec.getBaseQualities();
				
				
				if(cigar==null || cigar.isEmpty()) {
					SAMUtils.makeReadUnmapped(rec);
					sfw.addAlignment(rec);
					continue;
				}
				
				int readpos=0;
				for(final CigarElement ce:cigar) {
					final CigarOperator op = ce.getOperator();
					if(op.consumesReadBases()) {
						if(op.consumesReferenceBases()) {
							for(int x=0;x< ce.getLength();++x) {
								if(bases!=null) bases[readpos+x]= RESET_CHAR;
								if(quals!=null) quals[readpos+x]= RESET_QUAL;
								++nBasesMasked;
							}
						}
						readpos+=ce.getLength();
					}
				}
				++nRecordMasked;
				SAMUtils.makeReadUnmapped(rec);
				sfw.addAlignment(rec);
				}
			iter.close();
			sfw.close();
			progress.finish();
			LOG.info("done : reads masked "+nRecordMasked+"/"+nRecords+" Bases masked:"+nBasesMasked+"/"+nBases);
			return RETURN_OK;
			}
		catch(final Exception err)
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
		new SamMaskAlignedBases().instanceMainWithExit(args);
		}
	}
