package com.github.lindenb.jvarkit.tools.biostar;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;

@Deprecated /* use picard/RevertSam http://broadinstitute.github.io/picard/command-line-overview.html#RevertSam */
public class Biostar106668 extends AbstractCommandLineProgram
	{
	@Override
	public String getProgramDescription()
		{
		return "unmark duplicates. Deprecated: Use picard/RevertSam http://broadinstitute.github.io/picard/command-line-overview.html#RevertSam";
		}
	
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/Biostar106668";
		}

	@Override
	public void printOptions(PrintStream out)
		{
		out.print(" -b generates binary bam output "); 
		super.printOptions(out);
		}

	@Override
	public int doWork(String[] args)
		{
		boolean binary=false;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "b"))!=-1)
			{
			switch(c)
				{
				case 'b': binary= true;break;				
				default: 
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default: break;
						}
					}
				}
			}
		
		SamReader samReader=null;
		SAMRecordIterator iter=null;
		SAMFileWriter samWriter=null;
		long nConvert=0;
		try
			{
			SamFileReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
			if(opt.getOptInd()==args.length)
				{
				info("Reading sfomr stdin");
				samReader=SamFileReaderFactory.mewInstance().openStdin();
				}
			else if(opt.getOptInd()+1==args.length)
				{
				File filename=new File(args[opt.getOptInd()]);
				info("Reading from "+filename);
				samReader=SamFileReaderFactory.mewInstance().open(filename);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			SAMFileHeader header=samReader.getFileHeader();
			SAMProgramRecord prg=header.createProgramRecord();
			prg.setCommandLine(this.getProgramCommandLine());
			prg.setProgramName(this.getProgramName());
			prg.setProgramVersion(this.getVersion());
						
			SAMFileWriterFactory sfw=new SAMFileWriterFactory();
			sfw.setCreateIndex(false);
			sfw.setCreateMd5File(false);
			if(binary)
				{
				samWriter=sfw.makeBAMWriter(header,true, System.out);
				}
			else
				{
				samWriter=sfw.makeSAMWriter(header,true, System.out);
				}
			
			iter=samReader.iterator();
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(samReader.getFileHeader().getSequenceDictionary());
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
				progress.watch(rec);
				if(rec.getDuplicateReadFlag())
					{
					rec.setDuplicateReadFlag(false);
					++nConvert;
					}
				samWriter.addAlignment(rec);
				}
			progress.finish();
			info("Convert :"+nConvert);
			}
		catch (Exception e)
			{
			error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(samReader);
			CloserUtil.close(samWriter);
			}
		return 0;
		}

	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException
		{
		new Biostar106668().instanceMainWithExit(args);
		}
		

	}
