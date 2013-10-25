package com.github.lindenb.jvarkit.tools.tview;

import java.io.File;
import java.io.PrintWriter;

import com.github.lindenb.jvarkit.util.picard.IntervalUtils;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.util.CloserUtil;


public class TViewCmd extends CommandLineProgram
	{
	private static final Log LOG = Log.getInstance(TView.class);
    @Usage(programVersion="1.0")
    public String USAGE = getStandardUsagePreamble() + " java equivalent of samtools tview. ";
    @Option(shortName= StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Genome Reference",optional=false)
    public File REF=null;

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="A BAM file to process.")
    public File INPUT=null;
    @Option(shortName= "L", doc="Region to observe: chrom:start-end",optional=false)
    public String REGION="";
    
    
    @Override
    public String getVersion() {
    	return "1.0";
    	}
    
    @Override
	protected int doWork()
		{
    	IndexedFastaSequenceFile ref=null;
		PrintWriter out=new PrintWriter(System.out);
		SAMFileReader samReader=null;
		try {
	        samReader=new SAMFileReader(INPUT);
	        samReader.setValidationStringency(super.VALIDATION_STRINGENCY);
	        
			Interval interval=IntervalUtils.parseOne(samReader.getFileHeader().getSequenceDictionary(),REGION);
			if(interval==null)
				{
				LOG.error("Bad interval "+interval);
				return -1;
				}
			
			ref=new IndexedFastaSequenceFile(REF);
	        
	  
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
