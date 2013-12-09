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


import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionayProgress;
import com.github.lindenb.jvarkit.util.picard.SamSequenceRecordTreeMap;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.SequenceUtil;

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

    
	private int chrX_index=-1;
	private int chrY_index=-1;

    
    
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
    		boolean ok_pe_alignment=true;
			this.increment(Category.TOTAL);
			
			if(rec.getReadPairedFlag())
				{
				this.increment(Category.PAIRED);
				}
			else
				{
				ok_pe_alignment=false;
				}
			
			if(rec.getReadFailsVendorQualityCheckFlag())
				{
				ok_pe_alignment=false;
				this.increment(Category.FAIL_VENDOR_QUALITY);
				}
			
			if(rec.getReadUnmappedFlag())
				{
				ok_pe_alignment=false;
				this.increment(Category.UNMAPPED);
				return;
				}
			/* mapped below ************************************/
			
			this.increment(Category.MAPPED);
			
			if(rec.getReadPairedFlag() && rec.getProperPairFlag())
				{
				this.increment(Category.PROPER_PAIR);
				if(rec.getMappingQuality()>=QUAL)
					{
					this.increment(Category.PROPER_PAIR_HMQ);
					}
				}
			else
				{
				ok_pe_alignment=false;	
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
			else
				{
				ok_pe_alignment=false;
				}
			
			if(rec.getDuplicateReadFlag())
				{
				this.increment(Category.DUPLICATE);
				ok_pe_alignment=false;
				}
			
			if(rec.getMappingQuality()<QUAL)
				{
				this.increment(Category.FAIL_MAPPING_QUALITY);
				ok_pe_alignment=false;
				}
			
			
			
			if(ok_pe_alignment)
				{
				this.increment(Category.OK_FOR_PE_CALLING);
				if(chrX_index==rec.getReferenceIndex())
					{
					this.increment(Category.X);
					}
				else if(chrY_index==rec.getReferenceIndex())
					{
					this.increment(Category.Y);
					}

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
		PROPER_PAIR_HMQ,
		PLUS_STRAND,
		MINUS_STRAND,
		PRIMARY_ALIGNMENT,
		FAIL_MAPPING_QUALITY,
		DUPLICATE,
		FAIL_VENDOR_QUALITY,
		OK_FOR_PE_CALLING,
		X,Y;
		};
	
	@Override
	protected int doWork()
		{
		SAMFileReader samFileReader=null;
		//IntervalTreeMap<Interval> intervals=null;
		SamSequenceRecordTreeMap<Boolean> intervals=null;
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
				
				SAMSequenceDictionary samSequenceDictionary=null;
				
					
				for(File f:IN)
					{
					Map<String,Histogram> sample2hist=new HashMap<String, BamStats01.Histogram>();
					samFileReader=null;
					LOG.info("opening "+f);
					samFileReader=new SAMFileReader(f);
					samFileReader.setValidationStringency(super.VALIDATION_STRINGENCY);
					
					SAMSequenceDictionary currDict=samFileReader.getFileHeader().getSequenceDictionary();

					if(samSequenceDictionary==null)
						{
						samSequenceDictionary=currDict;
						}
					
					if(BEDILE!=null )
						{
						if(!SequenceUtil.areSequenceDictionariesEqual(currDict, samSequenceDictionary))
							{
							samFileReader.close();
							throw new IOException("incompatible sequence dictionaries. ("+f+")");
							}
							
						
						if(intervals==null)
							{
							intervals=new SamSequenceRecordTreeMap<Boolean>(currDict);
							LOG.info("opening "+BEDILE);
							Pattern tab=Pattern.compile("[\t]");
							String line;
							BufferedReader bedIn=IOUtils.openFileForBufferedReading(BEDILE);
							while((line=bedIn.readLine())!=null)
								{
								if(line.isEmpty() || line.startsWith("#")) continue;
								String tokens[]=tab.split(line,5);
								if(tokens.length<3) throw new IOException("bad bed line in "+line+" "+this.BEDILE);
								int seqIndex=currDict.getSequenceIndex(tokens[0]);
								if(seqIndex==-1)
									{
									throw new IOException("unknown chromosome from dict in  in "+line+" "+this.BEDILE);
									}
								int chromStart1= 1+Integer.parseInt(tokens[1]);//add one
								int chromEnd1= Integer.parseInt(tokens[2]);
								intervals.put(seqIndex, chromStart1, chromEnd1,Boolean.TRUE);
								}
							bedIn.close();
							LOG.info("done reading "+BEDILE);
							}
						}
					this.chrX_index=-1;
					this.chrY_index=-1;
					
					
					for(SAMSequenceRecord rec:currDict.getSequences())
						{
						String chromName=rec.getSequenceName().toLowerCase();
						if(chromName.equals("x") || chromName.equals("chrx"))
							{
							this.chrX_index=rec.getSequenceIndex();
							}
						else if(chromName.equals("y") || chromName.equals("chry"))
							{
							this.chrY_index=rec.getSequenceIndex();
							}
						}
					
					
					SAMSequenceDictionayProgress progess=new SAMSequenceDictionayProgress(currDict);
					SAMRecordIterator iter=samFileReader.iterator();
					while(iter.hasNext())
						{
						String sampleName=null;
						SAMRecord rec=iter.next();
						
						progess.watch(rec);
						
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
						if(rec.getReadUnmappedFlag())
							{
							continue;
							}
						
					
						
						if(!intervals.containsOverlapping(
									rec.getReferenceIndex(),
									rec.getAlignmentStart(),
									rec.getAlignmentEnd()
									))
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
