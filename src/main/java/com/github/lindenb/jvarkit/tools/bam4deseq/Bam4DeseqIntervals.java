package com.github.lindenb.jvarkit.tools.bam4deseq;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.github.lindenb.jvarkit.util.picard.IOUtils;


import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.SequenceUtil;

public class Bam4DeseqIntervals extends CommandLineProgram
	{
	private static final Log LOG=Log.getInstance(Bam4DeseqIntervals.class);
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" create a table for DESEQ with the number of reads within a sliding window for multiple BAMS.";
    
	@Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="BAM file to process",optional=false,minElements=1)
	public List<File> IN=new ArrayList<File>();
	@Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="output filename. Default stdout. ",optional=true)
	public File OUT=null;
	@Option(shortName="WIN",doc="size of the observed window. ",optional=true)
	public int WINDOW_SIZE=500;
	@Option(shortName="SHIFT",doc="shift window by SHIFT pb ",optional=true)
	public int WINDOW_SHIFT=250;
	@Option(shortName="COV",doc="ignore regions with NO coverage",optional=true)
	public boolean ONLY_COVERED=false;
	@Option(shortName="HEAD",doc="print header",optional=true)
	public boolean HEADER=true;
	

	@Override
	public String getVersion() {
		return "1.0";
		}
	
	@Override
	protected int doWork()
		{
		List<SAMFileReader> samFileReaders=new ArrayList<SAMFileReader>(IN.size());
		PrintWriter out=new PrintWriter(System.out);
		try {
			
			SAMSequenceDictionary ssDict=null;
			for(File in:IN)
				{
				LOG.info("opening "+in);
				SAMFileReader sfr=new SAMFileReader(in);
				sfr.setValidationStringency(super.VALIDATION_STRINGENCY);
				samFileReaders.add(sfr);
				
				SAMSequenceDictionary dict=sfr.getFileHeader().getSequenceDictionary();
				
				if(dict==null || dict.isEmpty())
					{
					LOG.error("No dictionary in "+in);
					return -1;
					}
				if(ssDict==null)
					{
					ssDict=dict;
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(dict, ssDict))
					{
					LOG.error("Not the same sequence dictionaries "+IN.get(0)+" vs "+in);
					return -1;
					}
				}
			if(OUT!=null)
				{
				out=new PrintWriter(IOUtils.openFileForBufferedWriting(OUT));
				}
			if(HEADER)
				{
				out.println("NAME");
				for(File in:IN)
					{
					out.print('\t');
					out.print(in);
					}
				out.println();
				}
			for(SAMSequenceRecord ssr:ssDict.getSequences())
				{
				LOG.info("scanning "+ssr.getSequenceName());
				int counts[]=new int[IN.size()];
				for(int i=1;i+this.WINDOW_SIZE <=ssr.getSequenceLength();i+=this.WINDOW_SHIFT)
					{
					boolean found=false;
					Arrays.fill(counts, 0);
					
					for(int t=0;t< samFileReaders.size();++t)
						{
						int count=0;
						SAMFileReader sfr=samFileReaders.get(t);
						SAMRecordIterator iter=sfr.queryOverlapping(ssr.getSequenceName(), i, i+this.WINDOW_SIZE);
						
						while(iter.hasNext())
							{
							SAMRecord rec=iter.next();
							if(rec.getDuplicateReadFlag()) continue;
							if(rec.getNotPrimaryAlignmentFlag()) continue;
							if(rec.getReadFailsVendorQualityCheckFlag()) continue;
							++count;
							}
						iter.close();
						if(count>0) found=true;
						counts[t]=count;
						}
					if(this.ONLY_COVERED && !found) continue;
					out.print(ssr.getSequenceName()+"_"+i+"_"+(i+this.WINDOW_SIZE));
					for(int count: counts)
						{
						out.print('\t');
						out.print(count);
						}
					out.println();
					}
				}
			return 0;
			} 
		catch (Exception e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			for(SAMFileReader sf:samFileReaders) sf.close();
			if(out!=null) out.flush();
			if(out!=null) out.close();
			}
		
		}
	
	public static void main(String[] args) {
		new Bam4DeseqIntervals().instanceMainWithExit(args);
		}
	}
