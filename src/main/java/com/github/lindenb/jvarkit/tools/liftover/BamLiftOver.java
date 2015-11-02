package com.github.lindenb.jvarkit.tools.liftover;

import java.util.Collection;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;

public class BamLiftOver extends AbstractBamLiftOver
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(BamLiftOver.class);

	@Override
	public Command createCommand() {
		return new MyCommand();
		}

	static private class MyCommand extends AbstractBamLiftOver.AbstractBamLiftOverCommand
		{    
		@Override
		protected Collection<Throwable> call(String inputName) throws Exception
			{
			SAMSequenceDictionary newDict=null;
			if(liftOverFile==null)
				{
				return wrapException("LiftOver file is undefined.");
				}
			if(REF==null)
				{
				return wrapException("REF undefined.");
				}
			
			SAMRecordIterator iter=null;
			SamReader sfr=null;
			SAMFileWriter sfw=null;
			try
				{
				newDict = new SAMSequenceDictionaryFactory().load(REF);
				sfr=openSamReader(inputName);
					
				
				LOG.info("Reading "+liftOverFile);
				LiftOver liftOver=new LiftOver(liftOverFile);
				liftOver.setLiftOverMinMatch(minMatch);
	
				SAMFileHeader headerIn=sfr.getFileHeader();
				
				SAMFileHeader headerOut=headerIn.clone();
				headerOut.setSortOrder(SortOrder.unsorted);
				
				SAMFileWriterFactory  sfwf=new SAMFileWriterFactory();
				sfwf.setCreateIndex(false);
				sfwf.setCreateMd5File(false);
				
				sfw = openSAMFileWriter(headerOut, true);
				
				iter=sfr.iterator();
				while(iter.hasNext())
					{
					SAMRecord rec=iter.next();
					SAMRecord copy=(SAMRecord)rec.clone();
					copy.setHeader(headerOut);
					StringBuilder sb=new StringBuilder();
					if(!rec.getReadUnmappedFlag())
						{
						String chrom=rec.getReferenceName();
						int pos=rec.getAlignmentStart();
						Interval interval=liftOver.liftOver(new Interval(chrom, pos,pos,rec.getReadNegativeStrandFlag(),null));
						if(interval!=null)
							{
							sb.append(chrom+":"+pos+":"+(rec.getReadNegativeStrandFlag()?"-":"+"));
							SAMSequenceRecord ssr=newDict.getSequence(interval.getContig());
							if(ssr==null)
								{
								sfr.close();
								sfr=null;
								sfw.close();
								return wrapException("the chromosome "+interval.getContig()+" is undefined in the sequence dict.");
								}
							copy.setReferenceName(ssr.getSequenceName());
							copy.setReferenceIndex(ssr.getSequenceIndex());
							copy.setAlignmentStart(interval.getStart());
							copy.setReadNegativeStrandFlag(interval.isNegativeStrand());
							
							}
						else
							{
							sb.append(".");
							copy.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
							copy.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
							copy.setAlignmentStart(0);
							copy.setMappingQuality(0);
							copy.setReadUnmappedFlag(true);
							copy.setSupplementaryAlignmentFlag(false);
							if(rec.getReadPairedFlag())
								{
								copy.setProperPairFlag(false);
								copy.setInferredInsertSize(0);
								}
							}
						}
					
					
					if(rec.getReadPairedFlag() && !rec.getMateUnmappedFlag())
						{
						sb.append("/");
						String chrom=rec.getMateReferenceName();
						int pos=rec.getMateAlignmentStart();
						Interval interval=liftOver.liftOver(new Interval(chrom, pos,pos,rec.getMateNegativeStrandFlag(),null));
						if(interval!=null)
							{
							sb.append(chrom+":"+pos+":"+(rec.getMateNegativeStrandFlag()?"-":"+"));
							SAMSequenceRecord ssr=newDict.getSequence(interval.getContig());
							if(ssr==null)
								{
								sfr.close();
								sfr=null;
								throw new SAMException("the chromosome "+interval.getContig()+" is undefined in the sequence dict.");
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
							copy.setMateReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
							copy.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
							copy.setMateAlignmentStart(0);
							copy.setMateUnmappedFlag(true);
							copy.setProperPairFlag(false);
							copy.setInferredInsertSize(0);
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
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new BamLiftOver().instanceMainWithExit(args);
		}

	}
