package com.github.lindenb.jvarkit.tools.liftover;

import java.io.File;
import java.io.PrintStream;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;

public class BamLiftOver extends AbstractCommandLineProgram
	{

	@Override
	public String getProgramDescription() {
		return "Lift-over a VCF file.";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfLiftOver";
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -f (chain-file) LiftOver file. Required.");
		out.println(" -m (double) lift over min-match. default:"+LiftOver.DEFAULT_LIFTOVER_MINMATCH);
		out.println(" -D (reference) indexed REFerence file for the new sequence dictionary. Required.");
		out.println(" -o (file.bam) file-out. Optional. Default stdout.");
		super.printOptions(out);
		}


	@Override
	public int doWork(String[] args)
		{
		
		File fileout=null;
		SAMSequenceDictionary newDict=null;
		double minMatch=LiftOver.DEFAULT_LIFTOVER_MINMATCH;
		File liftOverFile=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "f:m:X:D:o:"))!=-1)
			{
			switch(c)
				{
				case 'o': fileout=new File(opt.getOptArg());break;
				case 'D': 
					{
					try
						{
						newDict=new SAMSequenceDictionaryFactory().load(new File(opt.getOptArg()));
						}
					catch(Exception err)
						{
						error(err);
						return -1;
						}
					break;
					}
				case 'f': liftOverFile=new File(opt.getOptArg()); break;
				case 'm': minMatch=Double.parseDouble(opt.getOptArg()); break;
				default: 
					{
					switch(handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		if(liftOverFile==null)
			{
			error("LiftOver file is undefined.");
			return -1;
			}
		if(newDict==null)
			{
			error("New Sequence Dictionary file is undefined.");
			return -1;
			}
		SAMRecordIterator iter=null;
		SamReader sfr=null;
		SAMFileWriter sfw=null;
		try
			{
			if(opt.getOptInd()==args.length)
				{
				sfr=SamFileReaderFactory.mewInstance().openStdin();
				}
			else if(opt.getOptInd()+1==args.length)
				{
				File filein=new File(args[opt.getOptInd()]);
				sfr=SamFileReaderFactory.mewInstance().open(filein);
				}
			else
				{
				error("Illegal number of arguments");
				return -1;
				}
			info("Reading "+liftOverFile);
			LiftOver liftOver=new LiftOver(liftOverFile);
			liftOver.setLiftOverMinMatch(minMatch);

			SAMFileHeader headerIn=sfr.getFileHeader();
			
			SAMFileHeader headerOut=headerIn.clone();
			SAMProgramRecord prg=headerOut.createProgramRecord();
			prg.setProgramName(getProgramName());
			prg.setPreviousProgramGroupId(getProgramDescription());
			prg.setCommandLine(getProgramCommandLine());
			prg.setProgramName(getVersion());
			headerOut.setSortOrder(SortOrder.unsorted);
			
			SAMFileWriterFactory  sfwf=new SAMFileWriterFactory();
			sfwf.setCreateIndex(false);
			sfwf.setCreateMd5File(false);
			
			if(fileout==null)
				{
				info("Writing to stdout");
				sfwf.setTempDirectory(getTmpDirectories().get(0));
				sfw=sfwf.makeSAMWriter(headerOut, false, System.out);
				}
			else
				{
				info("Writing to "+fileout);
				super.addTmpDirectory(fileout);
				sfwf.setTempDirectory(super.getTmpDirectories().get(0));
				sfw=sfwf.makeSAMWriter(headerOut, false, fileout);
				}
			
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
						SAMSequenceRecord ssr=newDict.getSequence(interval.getSequence());
						if(ssr==null)
							{
							sfr.close();
							sfr=null;
							throw new SAMException("the chromosome "+interval.getSequence()+" is undefined in the sequence dict.");
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
						SAMSequenceRecord ssr=newDict.getSequence(interval.getSequence());
						if(ssr==null)
							{
							sfr.close();
							sfr=null;
							throw new SAMException("the chromosome "+interval.getSequence()+" is undefined in the sequence dict.");
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
			return 0;
			}
		catch(Exception err)
			{
			error(err);
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
