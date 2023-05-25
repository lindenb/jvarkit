/*
The MIT License (MIT)

Copyright (c) 2023 Pierre Lindenbaum

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
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassBamLauncher;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;

@Program(name="bamliftover",
	description="Lift-over a BAM file.",
	keywords={"bam","liftover"},
	generate_doc = true,
	jvarkit_amalgamion = true
	)
public class BamLiftOver extends OnePassBamLauncher
	{
	private static final Logger LOG = Logger.build(BamLiftOver.class).make();
	@Parameter(names={"-R2","--destination-dict"},description=DICTIONARY_SOURCE,required = true)
	private Path destDictionaryPath = null;


	@Parameter(names={"-f","--chain"},description="LiftOver file.",required = true)
	private File liftOverFile = null;

	@Parameter(names={"-m","--minmatch"},description="lift over min-match.")
	private double userMinMatch = LiftOver.DEFAULT_LIFTOVER_MINMATCH ;

	private LiftOver liftOver= null;
	private SAMFileHeader headerOut = null;
	private SAMSequenceDictionary newDict = null;
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int beforeSam() {
		final double minMatch=(this.userMinMatch<=0.0?LiftOver.DEFAULT_LIFTOVER_MINMATCH:this.userMinMatch);
		this.liftOver =new LiftOver(this.liftOverFile);
		this.liftOver.setLiftOverMinMatch(minMatch);
		this.newDict = SequenceDictionaryUtils.extractRequired(this.destDictionaryPath);
		return super.beforeSam();
		}
	
	@Override
	protected SAMFileHeader createOutputHeader(final SAMFileHeader headerIn) {
		this.headerOut =  super.createOutputHeader(headerIn);
		this.headerOut.setSortOrder(SortOrder.unsorted);
		return this.headerOut;
		}
	
	@Override
	protected Function<SAMRecord, List<SAMRecord>> createSAMRecordFunction() {
		return R->liftRead(R);
		}
	
	private List<SAMRecord> liftRead(final SAMRecord rec) {
		final SAMRecord copy;
		try {
			copy =(SAMRecord)rec.clone();
			}
		catch(final CloneNotSupportedException err) {
			throw new RuntimeException(err); 
			}
		copy.setHeader(this.headerOut);
		final StringBuilder sb=new StringBuilder();
		if(!rec.getReadUnmappedFlag())
			{
			final String chrom = rec.getReferenceName();
			final int pos=rec.getAlignmentStart();
			final Interval interval= this.liftOver.liftOver(new Interval(chrom, pos,pos,rec.getReadNegativeStrandFlag(),null));
			if(interval!=null)
				{
				sb.append(chrom+":"+pos+":"+(rec.getReadNegativeStrandFlag()?"-":"+"));
				final SAMSequenceRecord ssr= this.newDict.getSequence(interval.getContig());
				if(ssr==null)
					{
					throw new JvarkitException.ContigNotFoundInDictionary(interval.getContig(),newDict);
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
			final String chrom=rec.getMateReferenceName();
			final int pos=rec.getMateAlignmentStart();
			final Interval interval=liftOver.liftOver(new Interval(chrom, pos,pos,rec.getMateNegativeStrandFlag(),null));
			if(interval!=null)
				{
				sb.append(chrom+":"+pos+":"+(rec.getMateNegativeStrandFlag()?"-":"+"));
				final SAMSequenceRecord ssr=newDict.getSequence(interval.getContig());
				if(ssr==null)
					{
					throw new JvarkitException.ContigNotFoundInDictionary(interval.getContig(),newDict);
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
		rec.setAttribute(SAMTag.SA,null);
		if(sb.length()>0) copy.setAttribute("LO", sb.toString());
		return Collections.singletonList(copy);
		}
	

	public static void main(final String[] args)
		{
		new BamLiftOver().instanceMainWithExit(args);
		}

	}
