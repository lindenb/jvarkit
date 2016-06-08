
/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.bamstats01;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamSequenceRecordTreeMap;

import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;

public class BamStats01
	extends AbstractBamStats01
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(BamStats01.class);

	private PrintStream out=System.out;
	//private File bedFile=null;
	//private double minMappingQuality=30.0;
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
				if(rec.getMappingQuality()>=minMappingQuality &&
						rec.getMappingQuality()!=0 &&
						rec.getMappingQuality()!=255)
					{
					this.increment(Category.PROPER_PAIR_HMQ);
					}
				}
			else
				{
				ok_pe_alignment=false;	
				}		
			
			if(rec.getMappingQuality()==0 ||
				rec.getMappingQuality()==255)
				{
				this.increment(Category.BAD_MAPQ);
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
			
			if( !rec.isSecondaryOrSupplementary() && // bwa mem -P
				rec.getDuplicateReadFlag())
				{
				this.increment(Category.DUPLICATE);
				ok_pe_alignment=false;
				}
			
			if(rec.getMappingQuality()<minMappingQuality)
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
    	final Histogram2 histograms[]=new Histogram2[]{
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
		X,Y,
		BAD_MAPQ
		;
		};
	
	private BamStats01()
		{
		
		}
		
	private void run(final String filename,SamReader samFileReader) throws IOException
		{	
		Map<String,Histogram> sample2hist=new HashMap<String, BamStats01.Histogram>();
			
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
					LOG.info("opening "+this.bedFile);
					final Pattern tab=Pattern.compile("[\t]");
					String line;
					final BufferedReader bedIn=IOUtils.openFileForBufferedReading(bedFile);
					while((line=bedIn.readLine())!=null)
						{
						if(line.isEmpty() || line.startsWith("#")) continue;
						final String tokens[]=tab.split(line,5);
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
					LOG.info("done reading "+this.bedFile);
					}
				}
			this.chrX_index=-1;
			this.chrY_index=-1;
			
			
			for(final SAMSequenceRecord rec:currDict.getSequences())
				{
				final String chromName=rec.getSequenceName().toLowerCase();
				if(chromName.equals("x") || chromName.equals("chrx"))
					{
					this.chrX_index=rec.getSequenceIndex();
					}
				else if(chromName.equals("y") || chromName.equals("chry"))
					{
					this.chrY_index=rec.getSequenceIndex();
					}
				}
			
			
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(currDict);
			final SAMRecordIterator iter=samFileReader.iterator();
			while(iter.hasNext())
				{
				String sampleName=null;
				final SAMRecord rec=iter.next();
				
				progess.watch(rec);
				
				final SAMReadGroupRecord grp=rec.getReadGroup();
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
		
			for(final String sampleName: sample2hist.keySet())
				{
				final Histogram hist=sample2hist.get(sampleName);
				out.print(filename+"\t"+sampleName);
				
				for(final Category2 cat2: Category2.values())
					{
					for(final Category cat1: Category.values())//je je suis libertineuuh, je suis une cat1
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
	public Collection<Throwable> call() throws Exception {
		
		final List<String> args= new ArrayList<>(IOUtils.unrollFiles(this.getInputFiles())); 
		try {
			this.out = super.openFileOrStdoutAsPrintStream();
			
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
			
			final SamReaderFactory srf= super.createSamReaderFactory();
			
			if(args.isEmpty())
				{
				final SamReader r= srf.open(SamInputResource.of(stdin()));
				run("stdin",r);
				r.close();
				}
			else
				{
				for(final String filename:args)
					{
					LOG.info("Reading from "+filename);
					final SamReader sfr=srf.open(new File(filename));
					run(filename,sfr);
					sfr.close();
					}
				}
			out.flush();			
			this.out.flush();
			this.out.close();
			return RETURN_OK;
		} catch (Exception e) {
			return wrapException(e);
		} finally
			{
			CloserUtil.close(out);
			}
		}
	
	
	public static void main(String[] args)
		{
		new BamStats01().instanceMainWithExit(args);
		}
	}
