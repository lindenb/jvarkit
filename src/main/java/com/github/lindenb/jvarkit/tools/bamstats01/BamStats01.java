package com.github.lindenb.jvarkit.tools.bamstats01;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;


import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamSequenceRecordTreeMap;


import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileReader.ValidationStringency;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;

public class BamStats01
	extends AbstractCommandLineProgram
	{
	private PrintStream out=System.out;
	private File bedFile=null;
	private double QUAL=30.0;
	private int chrX_index=-1;
	private int chrY_index=-1;
	private SAMSequenceDictionary samSequenceDictionary=null;
	private SamSequenceRecordTreeMap<Boolean> intervals=null;
    
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
			
			/* count reads = unmapped+ primary align */
			if(rec.getReadUnmappedFlag() ||
				!rec.isSecondaryOrSupplementary())
				{
				this.increment(Category.UNMAPPED_PLUS_PRIMARY);
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
			
			if(rec.getSupplementaryAlignmentFlag())
				{
				this.increment(Category.SUPPLEMENTARY_ALIGNMENT);
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
		UNMAPPED_PLUS_PRIMARY,
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
		SUPPLEMENTARY_ALIGNMENT,//new in samSpec 2013-12
		X,Y;
		};
	
	private BamStats01()
		{
		
		}
		
	private void run(String filename,SAMFileReader samFileReader) throws IOException
		{	
		Map<String,Histogram> sample2hist=new HashMap<String, BamStats01.Histogram>();
			samFileReader.setValidationStringency(ValidationStringency.SILENT);
			
			SAMSequenceDictionary currDict=samFileReader.getFileHeader().getSequenceDictionary();

			if(samSequenceDictionary==null)
				{
				samSequenceDictionary=currDict;
				}
			
			if(this.bedFile!=null )
				{
				if(!SequenceUtil.areSequenceDictionariesEqual(currDict, samSequenceDictionary))
					{
					samFileReader.close();
					throw new IOException("incompatible sequence dictionaries."+filename);
					}
					
				
				if(intervals==null)
					{
					intervals=new SamSequenceRecordTreeMap<Boolean>(currDict);
					info("opening "+this.bedFile);
					Pattern tab=Pattern.compile("[\t]");
					String line;
					BufferedReader bedIn=IOUtils.openFileForBufferedReading(bedFile);
					while((line=bedIn.readLine())!=null)
						{
						if(line.isEmpty() || line.startsWith("#")) continue;
						String tokens[]=tab.split(line,5);
						if(tokens.length<3) throw new IOException("bad bed line in "+line+" "+this.bedFile);
						int seqIndex=currDict.getSequenceIndex(tokens[0]);
						if(seqIndex==-1)
							{
							throw new IOException("unknown chromosome from dict in  in "+line+" "+this.bedFile);
							}
						int chromStart1= 1+Integer.parseInt(tokens[1]);//add one
						int chromEnd1= Integer.parseInt(tokens[2]);
						intervals.put(seqIndex, chromStart1, chromEnd1,Boolean.TRUE);
						}
					bedIn.close();
					info("done reading "+this.bedFile);
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
			
			
			SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(currDict);
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
			progess.finish();
			samFileReader.close();
			samFileReader=null;
		
			for(String sampleName: sample2hist.keySet())
				{
				Histogram hist=sample2hist.get(sampleName);
				out.print(filename+"\t"+sampleName);
				
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
	
	@Override
	public String getProgramDescription()
		{
		return "Statistics about the reads in a BAM.";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/BamStats01";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println("-o (file) output. Default: stdout");
		out.println("-q (qual) min mapping quality: default: "+this.QUAL);
		out.println("-B (file) capture bed file.  Optional");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		File OUT=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:q:B:"))!=-1)
			{
			switch(c)
				{
				case 'o': OUT=new File(opt.getOptArg());break;
				case 'B': this.bedFile=new File(opt.getOptArg());break;
				case 'q': QUAL=Double.parseDouble(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		SAMFileReader samFileReader=null;		
		try
			{
			if(OUT!=null)
				{
				info("opening "+OUT+" for writing.");
				out=new PrintStream(IOUtils.openFileForWriting(OUT));
				}
			
			
			out.print("#Filename\tSample");
			for(Category2 cat2: Category2.values())
				{
				for(Category cat1: Category.values())
					{
					out.print("\t"+cat2+"_"+cat1);//je je suis libertineuuh, je suis une cat1
					}
				if(bedFile==null) break;
				}
			out.println();
			
			
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				SAMFileReader sfr=new SAMFileReader(System.in);
				run("stdin",sfr);
				sfr.close();
				}
			else
				{
				for(int i=opt.getOptInd();i< args.length;++i)
					{
					String filename=args[i];
					info("Reading from "+filename);
					SAMFileReader sfr=new SAMFileReader(new File(filename));
					run(filename,sfr);
					sfr.close();
					}
				}
			out.flush();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(samFileReader);
			CloserUtil.close(out);
			}
		}
	
	public static void main(String[] args)
		{
		new BamStats01().instanceMainWithExit(args);
		}
	}
