package com.github.lindenb.jvarkit.tools.bamstats01;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;


import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.IOUtils;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class BamStats01
	extends AbstractCommandLineProgram
	{
	private static final Log LOG=Log.getInstance(BamStats01.class);
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" Statistics about the reads in a BAM. ";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="BAM files to process.",
    		optional=false,
    		minElements=1)
	public List<File> IN=new ArrayList<File>();

    @Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME,
    		doc="Ouput. Default stdout. ",
    		optional=true)
	public File OUT=null;

    
    
    @Option(shortName="Q", doc="Default treshold quality", optional=true)
	public double QUAL=30.0;
	
    @Option(shortName= "BED", doc="BED File.",
    		optional=true)
	public File BEDILE=null;

    
    private  class Histogram2
    	{
    	long counts[];
    	
    	Histogram2()
    		{
    		counts=new long[Category.values().length];
    		Arrays.fill(counts, 0L);
    		}
    	void increment(Category cat)
    		{
    		counts[cat.ordinal()]++;
    		}
    	void watch(SAMRecord rec)
    		{
			this.increment(Category.TOTAL);
			
			if(rec.getReadPairedFlag())
				{
				this.increment(Category.PAIRED);
				
				}
			
			if(rec.getReadFailsVendorQualityCheckFlag())
				{
				this.increment(Category.FAIL_VENDOR_QUALITY);
				}
			
			if(rec.getReadUnmappedFlag())
				{
				this.increment(Category.UNMAPPED);
				return;
				}
			/* mapped below ************************************/
			
			this.increment(Category.MAPPED);
			
			if(rec.getReadPairedFlag() && rec.getProperPairFlag())
				{
				this.increment(Category.PROPER_PAIR);
				}
			
			
			
			if(rec.getReadNegativeStrandFlag())
				{
				this.increment(Category.MINUS_STRAND);
				}
			else
				{
				this.increment(Category.PLUS_STRAND);
				}
			
			
			if(!rec.getNotPrimaryAlignmentFlag())
				{
				this.increment(Category.PRIMARY_ALIGNMENT);
				}
			
			
			if(rec.getDuplicateReadFlag())
				{
				this.increment(Category.DUPLICATE);
				}
			
			if(rec.getMappingQuality()<QUAL)
				{
				this.increment(Category.FAIL_MAPPING_QUALITY);
				}

    		}
    	
    	}
    
    private  class Histogram
    	{
    	Histogram2 histograms[]=new Histogram2[]{
    			new Histogram2(),//ALL
    			new Histogram2(),//in capture
    			new Histogram2()//off capture
    			};
    	}

    
    private enum Category2
    	{
    	ALL,IN_TARGET,OFF_TARGET
    	}
	
	
	private enum Category
		{
		TOTAL,
		PAIRED,
		UNMAPPED,
		MAPPED,
		PROPER_PAIR,
		PLUS_STRAND,
		MINUS_STRAND,
		PRIMARY_ALIGNMENT,
		FAIL_MAPPING_QUALITY,
		DUPLICATE,
		FAIL_VENDOR_QUALITY;
		
		};
	
	@Override
	protected int doWork()
		{
		SAMFileReader samFileReader=null;
		IntervalTreeMap<Interval> intervals=null;
		
		PrintStream out=System.out;
		
		
		
			try
				{
				if(OUT!=null)
					{
					LOG.info("opening "+OUT+" for writing.");
					out=new PrintStream(IOUtils.openFileForWriting(OUT));
					}
				
				
				out.print("#Filename\tSample");
				for(Category2 cat2: Category2.values())
					{
					for(Category cat1: Category.values())
						{
						out.print("\t"+cat2+"_"+cat1);//je je suis libertineuuh, je suis une cat1
						}
					if(BEDILE==null) break;
					}
				out.println();
				
				
				if(BEDILE!=null)
					{
					intervals=new IntervalTreeMap<Interval>();
					LOG.info("opening "+BEDILE);
					Pattern tab=Pattern.compile("[\t]");
					String line;
					BufferedReader bedIn=IOUtils.openFileForBufferedReading(BEDILE);
					while((line=bedIn.readLine())!=null)
						{
						if(line.isEmpty() || line.startsWith("#")) continue;
						String tokens[]=tab.split(line,5);
						if(tokens.length<3) throw new IOException("bad bed line in "+line+" "+this.BEDILE);
						String chrom=tokens[0];
						int chromStart1= 1+Integer.parseInt(tokens[1]);//add one
						int chromEnd1= Integer.parseInt(tokens[2]);
						Interval interval=new Interval(chrom, chromStart1, chromEnd1);
						intervals.put(interval, interval);
						}
					bedIn.close();
					}
					
				for(File f:IN)
					{
					Map<String,Histogram> sample2hist=new HashMap<String, BamStats01.Histogram>();
					samFileReader=null;
					LOG.info("opening "+f);
					samFileReader=new SAMFileReader(f);
					samFileReader.setValidationStringency(super.VALIDATION_STRINGENCY);
					SAMRecordIterator iter=samFileReader.iterator();
					while(iter.hasNext())
						{
						String sampleName=null;
						SAMRecord rec=iter.next();
						SAMReadGroupRecord grp=rec.getReadGroup();
						if(grp!=null)
							{
							sampleName=grp.getSample();
							}
						
						if(sampleName==null || sampleName.isEmpty()) sampleName="undefined";
						
						Histogram hist=sample2hist.get(sampleName);
						if(hist==null)
							{
							hist=new Histogram();
							sample2hist.put(sampleName, hist);
							}
						
						hist.histograms[Category2.ALL.ordinal()].watch(rec);
						
						if(intervals==null) continue;
						if(rec.getReadUnmappedFlag()) continue;
						
						if(intervals.getOverlapping(new Interval(
									rec.getReferenceName(),
									rec.getAlignmentStart(),
									rec.getAlignmentEnd()
									)).isEmpty()
							)
							{
							hist.histograms[Category2.OFF_TARGET.ordinal()].watch(rec);
							}		
						else
							{
							hist.histograms[Category2.IN_TARGET.ordinal()].watch(rec);
							}
						}
					samFileReader.close();
					samFileReader=null;
					for(String sampleName: sample2hist.keySet())
						{
						Histogram hist=sample2hist.get(sampleName);
						out.print(f.getPath()+"\t"+sampleName);
						
						for(Category2 cat2: Category2.values())
							{
							for(Category cat1: Category.values())//je je suis libertineuuh, je suis une cat1
								{
								out.print("\t");
								out.print(hist.histograms[cat2.ordinal()].counts[cat1.ordinal()]);
								}
							if(intervals==null) break;
							}
						out.println();
						}
					}
				out.flush();
		        }
			catch(Exception err)
				{
				LOG.error(err, ""+err.getMessage());
				return -1;
				}
			finally
				{
				if(samFileReader!=null) samFileReader.close();
				out.flush();
				if(OUT!=null) {out.close();}
				}
		
		return 0;
		}
	
	public static void main(String[] args)
		{
		new BamStats01().instanceMainWithExit(args);
		}
	}
