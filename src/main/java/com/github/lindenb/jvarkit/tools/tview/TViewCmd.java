package com.github.lindenb.jvarkit.tools.tview;

import java.io.File;
import java.io.PrintWriter;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.IntervalUtils;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileReader.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;


public class TViewCmd extends AbstractCommandLineProgram
	{

    
    @Override
	public String getProgramDescription() {
		return "java equivalent of samtools tview";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println("-p (chrom:pos) region");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		String region=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args, getGetOptDefault()+"r:"))!=-1)
			{
			switch(c)
				{
				case 'r': region=opt.getOptArg();break;
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
		
		File faidx=null;
		File bamFile=null;
		
		if(opt.getOptInd()+1==args.length )
			{
			bamFile=new File(args[opt.getOptInd()]);
			}
		else if(opt.getOptInd()+2==args.length )
			{
			bamFile=new File(args[opt.getOptInd()]);
			faidx=new File(args[opt.getOptInd()+1]);
			}
		else
			{
			System.err.println("Illegal Number of arguments.");
			return -1;
			}
		
		
	
		
    	IndexedFastaSequenceFile ref=null;
		PrintWriter out=new PrintWriter(System.out);
		SAMFileReader samReader=null;
		try {
			info("opening "+bamFile);
	        samReader=new SAMFileReader(bamFile);
	        samReader.setValidationStringency(ValidationStringency.LENIENT);
	        
			Interval interval=IntervalUtils.parseOne(
					samReader.getFileHeader().getSequenceDictionary(),
					region
					);
			if(interval==null)
				{
				error("Bad interval "+interval);
				return -1;
				}
			if(faidx!=null)
				{
				ref=new IndexedFastaSequenceFile(faidx);
				}
	  
	        TViewHandler handler=new AsciiHandler();
	        new TView().execute(samReader, ref, interval, handler);
			} 
		catch (Exception e) {
			e.printStackTrace();
			return -1;
			}
		finally
			{
			out.flush();
			CloserUtil.close(samReader);
			CloserUtil.close(ref);
			}
		return 0;
		}

public static void main(String[] args)
	{
	new TViewCmd().instanceMainWithExit(args);
	}
}
