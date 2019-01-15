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
package com.github.lindenb.jvarkit.tools.liftover;
/**
BEGIN_DOC



END_DOC
 */
import java.io.File;
import java.util.List;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
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


import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

@Program(name="bamliftover",
	description="Lift-over a BAM file.",
	keywords={"bam","liftover"}
		)
public class BamLiftOver extends Launcher
	{
	private static final Logger LOG = Logger.build(BamLiftOver.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-f","--chain"},description="LiftOver file. Require")
	private File liftOverFile = null;

	@Parameter(names={"-m","--minmatch"},description="lift over min-match. default:-1 == use default value from htsjdk LiftOver.DEFAULT_LIFTOVER_MINMATCH")
	private double userMinMatch = -1 ;

	@Parameter(names={"-D","-R","--reference"},description="indexed REFerence file for the new sequence dictionary. Required")
	private File faidx = null;

	@ParametersDelegate
	private WritingBamArgs writingBamArgs =new WritingBamArgs();
	
	@Override
	public int doWork(final List<String> args) {
		final double minMatch=(this.userMinMatch<=0.0?LiftOver.DEFAULT_LIFTOVER_MINMATCH:this.userMinMatch);
		if(this.liftOverFile==null)
			{
			LOG.error("LiftOver file is undefined.");
			return -1;
			}
		if(this.faidx==null)
			{
			LOG.error("New Sequence Dictionary file is undefined.");
			return -1;
			}
		SAMRecordIterator iter=null;
		SamReader sfr=null;
		SAMFileWriter sfw=null;
		try
			{
			LOG.info("Reading "+liftOverFile);
			LiftOver liftOver=new LiftOver(liftOverFile);
			liftOver.setLiftOverMinMatch(minMatch);

			
			final SAMSequenceDictionary newDict=SAMSequenceDictionaryExtractor.extractDictionary(faidx);
			
			
			sfr=super.openSamReader(oneFileOrNull(args));
			

			final SAMFileHeader headerIn=sfr.getFileHeader();
			final SAMFileHeader headerOut=headerIn.clone();
			headerOut.setSortOrder(SortOrder.unsorted);
			
			sfw = this.writingBamArgs.openSAMFileWriter(outputFile,headerOut, true);
			
			
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
							LOG.error("the chromosome "+interval.getContig()+" is undefined in the sequence dict.");
							return -1;
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
							LOG.error("the chromosome "+interval.getContig()+" is undefined in the sequence dict.");
							return -1;
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
		new BamLiftOver().instanceMainWithExit(args);
		}

	}
