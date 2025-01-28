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
package com.github.lindenb.jvarkit.tools.liftover;
import java.io.IOException;
/**
BEGIN_DOC



END_DOC
 */
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassBamLauncher;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.ucsc.LiftOverChain;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMException;
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
	modificationDate = "20250114",
	generate_doc = true,
	jvarkit_amalgamion = true
	)
public class BamLiftOver extends OnePassBamLauncher
	{
	private static final Logger LOG = Logger.build(BamLiftOver.class).make();
	@Parameter(names={"-R2","--destination-dict"},description=DICTIONARY_SOURCE,required = true)
	private Path destDictionaryPath = null;


	@Parameter(names={"-f","--chain"},description=LiftOverChain.OPT_DESC,required = true)
	private String liftOverFile = null;

	@Parameter(names={"-m","--minmatch"},description="lift over min-match.")
	private double userMinMatch = LiftOver.DEFAULT_LIFTOVER_MINMATCH ;

	@Parameter(names={"--unmapped"},description="discard unmapped reads/unlifted reads")
	private boolean discard_unmapped = false ;

	@Parameter(names={"--save-position"},description="Save original position in SMA attribute")
	private boolean  save_original_position = false ;
	@Parameter(names={"--drop-seq"},description="drop SEQ and QUAL")
	private boolean  drop_seq_and_qual = false ;

	
	private LiftOver liftOver= null;
	private SAMFileHeader headerOut = null;
	private SAMSequenceDictionary newDict = null;
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int beforeSam() {
		this.newDict = SequenceDictionaryUtils.extractRequired(this.destDictionaryPath);
		if(StringUtils.isBlank(this.liftOverFile)) {
			LOG.error("chain is blank");
			return -1;
			}
		return super.beforeSam();
		}
	
	@Override
	protected SAMFileHeader createOutputHeader(final SAMFileHeader headerIn) {
		
		try {
			final SAMSequenceDictionary dictIn=SequenceDictionaryUtils.extractRequired(headerIn);
			final double minMatch=(this.userMinMatch<=0.0?LiftOver.DEFAULT_LIFTOVER_MINMATCH:this.userMinMatch);
			this.liftOver = LiftOverChain.load(this.liftOverFile,dictIn,this.newDict);
			this.liftOver.setLiftOverMinMatch(minMatch);
			}
		catch(IOException errr) {
			throw new SAMException(errr);
			}
	

		this.headerOut =  super.createOutputHeader(headerIn);
		this.headerOut.setSequenceDictionary(this.newDict);
		this.headerOut.setSortOrder(SortOrder.unsorted);
		return this.headerOut;
		}
	
	@Override
	protected Function<SAMRecord, List<SAMRecord>> createSAMRecordFunction() {
		return R->liftRead(R);
		}
	
	private SAMRecord dropSeqAndQual(final SAMRecord rec) {
		if(this.drop_seq_and_qual) {
			rec.setReadBases(SAMRecord.NULL_QUALS);
			rec.setBaseQualities(SAMRecord.NULL_QUALS);
			}
		return rec;
		}
	
	private void makeReadUnmapped(final SAMRecord rec) {
		SAMUtils.makeReadUnmapped(dropSeqAndQual(rec));
		
	}
	
	private List<SAMRecord> liftRead(final SAMRecord source) {
		final SAMRecord copy;
		try {
			copy =(SAMRecord)source.clone();
			}
		catch(final CloneNotSupportedException err) {
			throw new RuntimeException(err); 
			}
		copy.setHeader(this.headerOut);
		final StringBuilder sb =(save_original_position?new StringBuilder():null);
		if(!source.getReadUnmappedFlag())
			{
			final String chrom = source.getReferenceName();
			final int pos=source.getAlignmentStart();
			Interval interval= this.liftOver.liftOver(new Interval(chrom, pos,pos,source.getReadNegativeStrandFlag(),null));
			if(interval!=null && this.newDict.getSequence(interval.getContig())==null) interval=null;
			if(interval!=null)
				{
				if(sb!=null) sb.append(chrom+":"+pos+":"+(source.getReadNegativeStrandFlag()?"-":"+"));
				final SAMSequenceRecord ssr= this.newDict.getSequence(interval.getContig());
				if(ssr==null)
					{
					throw new JvarkitException.ContigNotFoundInDictionary(interval.getContig(),newDict);
					}
				copy.setReferenceName(ssr.getSequenceName());
				copy.setReferenceIndex(ssr.getSequenceIndex());
				copy.setAlignmentStart(interval.getStart());
				copy.setReadNegativeStrandFlag(interval.isNegativeStrand());
				
				if(drop_seq_and_qual) {
					copy.setReadBases(SAMRecord.NULL_QUALS);
					copy.setBaseQualities(SAMRecord.NULL_QUALS);
					}
				else if(source.getReadNegativeStrandFlag()!=copy.getReadNegativeStrandFlag()) {
					copy.setReadString(AcidNucleics.reverseComplement(source.getReadString()));
					
					final byte qual[]= source.getBaseQualities();
					final byte quals2[]=  new byte[qual.length];
					for(int i=0;i< qual.length;++i) {
						quals2[i]=qual[(qual.length-1)-i];
					}
					copy.setBaseQualities(quals2);
					}
				}
			else
				{
				if(discard_unmapped) return Collections.emptyList();
				if(sb!=null) sb.append(".");
				this.makeReadUnmapped(copy);
				}
			}
		else
			{
			if(discard_unmapped) return Collections.emptyList();
			this.makeReadUnmapped(copy);
			}
		
		
		if(source.getReadPairedFlag())
			{
			if(source.getMateReferenceIndex()!=SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) /* source.getMateUnmappedFlag() */{
				if(sb!=null)sb.append("/");
				final String chrom = source.getMateReferenceName();
				if(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(chrom)) throw new SAMException("illegal state "+source.getReadName());
				final int pos = source.getMateAlignmentStart();
				Interval interval=liftOver.liftOver(new Interval(chrom, pos,pos,source.getMateNegativeStrandFlag(),null));
				if(interval!=null && this.newDict.getSequence(interval.getContig())==null) {
					interval=null;
					}
				if(interval!=null)
					{
					if(sb!=null) sb.append(chrom+":"+pos+":"+(source.getMateNegativeStrandFlag()?"-":"+"));
					final SAMSequenceRecord ssr= this.newDict.getSequence(interval.getContig());
					if(ssr==null)
						{
						throw new JvarkitException.ContigNotFoundInDictionary(interval.getContig(),this.newDict);
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
					if(sb!=null) sb.append(".");
					copy.setMateUnmappedFlag(true);
					copy.setMateReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
					copy.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
					copy.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
					copy.setProperPairFlag(false);
					copy.setInferredInsertSize(0);
					}
				}
			else //mate unmapped BUT mate-reference is defined
				{
				if(copy.getReadUnmappedFlag()) {
					copy.setMateReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
					copy.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
					copy.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
					copy.setProperPairFlag(false);
					copy.setInferredInsertSize(0);
					}
				else
					{
					copy.setMateReferenceName(copy.getReferenceName());
					copy.setMateReferenceIndex(copy.getReferenceIndex());
					copy.setMateAlignmentStart(copy.getAlignmentStart());
					copy.setProperPairFlag(false);
					copy.setInferredInsertSize(0);
					}
				}
			}
		/*
		if(
			(copy.getReferenceName()!=null && !copy.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) && this.newDict.getSequence(copy.getReferenceName())==null) ||
			(copy.getMateReferenceName()!=null && !copy.getMateReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) && this.newDict.getSequence(copy.getMateReferenceName())==null)
			) {
			throw new IllegalStateException("\nsource: "+source.getSAMString()+"\ndest "+copy.getSAMString()+"\n");
		}*/
		
		copy.setAttribute(SAMTag.SA,null);
		if(sb!=null && sb.length()>0) copy.setAttribute("LO", sb.toString());
		dropSeqAndQual(copy);
		return Collections.singletonList(copy);
		}
	

	public static void main(final String[] args)
		{
		new BamLiftOver().instanceMainWithExit(args);
		}

	}
