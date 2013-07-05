package com.github.lindenb.jvarkit.tools.bamstats01;

import java.io.File;
import java.util.Map;
import java.util.logging.Logger;

import net.sf.picard.util.IntervalTree;
import com.github.lindenb.jvarkit.tools.cmpbams.CompareBams2;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.metrics.MetricBase;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.util.Histogram;
import net.sf.picard.util.RExecutor;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class BamStats01 extends CommandLineProgram
	{
	private static final Logger LOG=Logger.getLogger(CompareBams2.class.getName());
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
		try
			{
			SAMRecordIterator iter=samFileReader.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
				hist.increment(Category.ALL);
				if(rec.getReadUnmappedFlag())
					{
					hist.increment(Category.UNMAPPED);
					continue;
					}
				hist.increment(Category.MAPPED);
				if(rec.getDuplicateReadFlag())
					{
					hist.increment(Category.DUPLICATE);
					}
				if(rec.getReadFailsVendorQualityCheckFlag())
					{
					hist.increment(Category.FAIL_VENDOR_QUALITY);
					}
				if(QUAL!=null && rec.getMappingQuality()<QUAL)
					{
					hist.increment(Category.FAIL_MAPPING_QUALITY);
					}
				net.sf.picard.util.IntervalTree<Boolean> iTree=null;
				
				if(bed!=null &&
					((iTree=bed.get(rec.getReferenceIndex()))!=null) &&
					iTree.overlappers(rec.getAlignmentStart(), rec.getAlignmentEnd()).hasNext())
					{
					hist.increment(Category.IN_CAPTURE);
					}
				}
			final MetricsFile<Metrics,Category> metrics = getMetricsFile();
	        metrics.addHistogram(hist);
	        metrics.write(OUT);

	        if (hist.isEmpty() && hist.isEmpty()) {
	        	LOG.warning("No valid bases found in input file. No plot will be produced.");
	        }
	        else {
	            // Now run R to generate a chart
	            final int rResult = RExecutor.executeFromClasspath(
	                    "net/sf/picard/analysis/qualityScoreDistribution.R",
	                    OUT.getAbsolutePath(),
	                    OUT.getAbsolutePath(),
	                    IN.getName());

	            if (rResult != 0) {
	                throw new PicardException("R script qualityScoreDistribution.R failed with return code " + rResult);
	            }
	        
		
	        	}
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
