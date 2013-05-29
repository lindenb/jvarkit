package com.github.lindenb.jvarkit.tools.sam4weblogo;

import java.io.File;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;


import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.picard.util.Log;
import net.sf.picard.util.SamLocusIterator;
import net.sf.picard.util.SamLocusIterator.RecordAndOffset;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

public class SAM4WebLogo extends CommandLineProgram
	{
	@SuppressWarnings("unused")
	private static final Log log = Log.getInstance(SAM4WebLogo.class);
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Sequence logo for different alleles or generated from SAM/BAM http://www.biostars.org/p/73021";
	
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="A BAM file to process.")
    public File INPUT=null;
    @Option(shortName= "R", doc="Region to observe: chrom:start-end",optional=false)
    public String REGION="";
    
    
    
    @Override
	protected int doWork()
		{
		PrintWriter out=new PrintWriter(System.out);
		SAMFileReader samReader=null;
		SamLocusIterator slit=null;
		Iterator<SamLocusIterator.LocusInfo> iter=null;
		try {
			Interval interval=parseInterval(REGION);
			if(interval==null) return -1;
	        samReader=new SAMFileReader(INPUT);
	        samReader.setValidationStringency(super.VALIDATION_STRINGENCY);
	        Map<SAMRecord, StringBuilder> record2seq =new LinkedHashMap<SAMRecord, StringBuilder>();
            IntervalList  iL=new  IntervalList(samReader.getFileHeader());
            iL.add(interval);
	        slit=new  SamLocusIterator(samReader,iL,true);
	        slit.setEmitUncoveredLoci(true);
	        int length=0;
	        for(iter=slit.iterator();
                    iter.hasNext();
                    )
                {
                SamLocusIterator.LocusInfo  locusInfo=iter.next();
                if(locusInfo.getPosition() < interval.getStart() ) continue;
                if(locusInfo.getPosition() > interval.getEnd() ) continue;
                for(RecordAndOffset rao: locusInfo.getRecordAndPositions())
                	{
                	SAMRecord rec=rao.getRecord();
                	StringBuilder b=record2seq.get(rec);
                	if(b==null)
                		{
                		b=new StringBuilder();
                		while(b.length()<length) b.append("-");
                		record2seq.put(rec,b);
                		}
                	char c=(char)rao.getReadBase();
                	b.append(c);
                	}
                ++length;
                for(SAMRecord rec: record2seq.keySet())
                	{
                	StringBuilder b=record2seq.get(rec);
                	while(b.length()<length) b.append("-");
                	}
               
                }
	        for(SAMRecord rec: record2seq.keySet())
	        	{
	        	StringBuilder b=record2seq.get(rec);
	        	out.print(">"+rec.getReadName());
	        	if(rec.getReadPairedFlag())
	        		{
	        		if(rec.getFirstOfPairFlag()) out.print("/1");
	        		if(rec.getSecondOfPairFlag()) out.print("/2");
	        		}
	        	out.println();
	        	out.println(b);
	        	}
	        
			} 
		catch (Exception e) {
			e.printStackTrace();
			return -1;
			}
		finally
			{
			if(slit!=null) slit.close();
			if(samReader!=null) samReader.close();
			out.flush();
			}
		return 0;
		}

private static Interval parseInterval(String reg)
	{
	try
			{	 
				int colon = reg.indexOf(':');
				if(colon<1) throw new IllegalArgumentException("bad region "+reg);
				int hyphen = reg.indexOf('-');
		
				String s=reg.substring(0,colon);
				int start= Integer.parseInt(reg.substring(colon+1,hyphen));
				int end=Integer.parseInt(reg.substring(hyphen+1));
				return new Interval(s, start, end);
			}
			catch(Exception err)
			{
				System.err.println("bad interval "+reg);
				return null;
			}
		}
public static void main(String[] args)
	{
	new SAM4WebLogo().instanceMainWithExit(args);
	}
}
