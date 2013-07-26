package com.github.lindenb.jvarkit.tools.bamstats01;

import java.io.File;
import java.util.Map;
import java.util.logging.Logger;

import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.metrics.MetricBase;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.util.Histogram;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecordIterator;

public class BamStats01 extends AbstractCommandLineProgram
	{
	private static final Logger LOG=Logger.getLogger(BamStats01.class.getName());
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"Compare two or more BAM files.";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="BAM files to process.",
    		optional=false)
	public File IN=null;

    @Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="BAM files to process.",
    		optional=false)
	public File OUT=null;

	public Double QUAL=null;
	
	public static class Metrics extends MetricBase
		{
	
		}
	
	public static class BamStats01Metrics
		extends MetricsFile<Metrics, Category>
		{
	
		}
	
	private enum Category {
		ALL,
		UNMAPPED,
		MAPPED,
		FAIL_MAPPING_QUALITY,
		DUPLICATE,
		FAIL_VENDOR_QUALITY,
		IN_CAPTURE
		};
	@Override
	protected int doWork()
		{
		Map<Integer,net.sf.picard.util.IntervalTree<Boolean>> bed=null;
		Histogram<Category> hist=new Histogram<BamStats01.Category>();
		SAMFileReader samFileReader=null;
		
		BamStats01Report report=new BamStats01Report(samFileReader.getFileHeader());
		try
			{
			SAMRecordIterator iter=samFileReader.iterator();
			report.addAlignment(iter.next());
	        }
		catch(Exception err)
			{
			return -1;
			}
		finally
			{
			
			}
		
		return 0;
		}
	
	public static void main(String[] args)
		{
		new BamStats01().instanceMainWithExit(args);
		}
	}
