
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
* 2016 : creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.util.List;

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
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**

BEGIN_DOC


Supplementary and secondary reads are ignored.
Output Bam is unsorted and reads are all unmapped.



END_DOC
*/
@Program(name="sammaskalignedbases",description="Mask bases aligned on Reference.")
public class SamMaskAlignedBases
	extends Launcher
	{
	private static final Logger LOG = Logger.build(SamMaskAlignedBases.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@ParametersDelegate
	private WritingBamArgs writingBamArgs =new WritingBamArgs();
	
	
	public SamMaskAlignedBases()
	{
	}
	
	@Override
	public int doWork(List<String> args) {
	
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
			sfr = openSamReader(oneFileOrNull(args));
			
			final SAMFileHeader header1=sfr.getFileHeader();
			
			if(header1==null)
				{
				LOG.error("File header missing");
				return -1;
				}
			final SAMFileHeader header2=header1.clone();
			header2.addComment(getProgramName()+":"+getVersion()+":"+getProgramCommandLine());
			header2.setSortOrder(SortOrder.unsorted);
			
			sfw =  this.writingBamArgs.openSAMFileWriter(outputFile,header2, true);
			
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
			LOG.error(err);
			return -1;
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
