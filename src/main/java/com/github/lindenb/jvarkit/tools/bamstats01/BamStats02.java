
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

import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.illumina.ShortReadName;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;

public class BamStats02
	extends AbstractBamStats02
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(BamStats02.class);
	private IntervalTreeMap<Boolean> intervals=null;
    
   
	private static boolean isEmpty(String s)
		{
		return s==null || s.trim().isEmpty() || s.equals(".");
		}
	
    /* visible from BamStats02View */ static enum STRING_PROPS {
    	filename,samplename,
    	chromosome,mate_chromosome,
    	platform,platformUnit,
    	instrument,flowcell,
    	library
    	};
    /* visible from BamStats02View */ static enum INT_PROPS {lane,run,mapq,inTarget};
    
    /** tuple of properties */
    /* visible from BamStats02View */ static class Category
    	{
    	String _strings[]=new String[STRING_PROPS.values().length];
    	int _ints[]=new int[INT_PROPS.values().length];
    	int flag=-1;
    	
    	Category()
    		{	
    		Arrays.fill(_strings, ".");
    		Arrays.fill(_ints, -1);
    		}
    	
    	private void print(PrintWriter pw)
    		{
    		boolean first=true;
    		for(String s:_strings)
    			{
    			if(!first) pw.print("\t");
    			first=false;
    			pw.print(s);
    			}
    		for(int s:_ints)
				{
    			pw.print("\t");
				pw.print(s);
				}
    		
    		for(final SAMFlag flg:SAMFlag.values())
    			{
    			pw.print("\t");
    			pw.print(flg.isSet(this.flag)?1:0);
    			}
    		}
    	
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result
					+ Arrays.hashCode(_strings);
			result = prime * result
					+ Arrays.hashCode(_ints);
			result = prime * result + flag;
			return result;
			}
		
		
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Category other = (Category) obj;
			
			if (flag != other.flag)
			if(!Arrays.equals(this._ints, other._ints)) return false;
			if(!Arrays.equals(this._strings, other._strings)) return false;			
			return true;
			}
    	public String getChrom()
    		{
    		String s= _strings[STRING_PROPS.chromosome.ordinal()];
    		return isEmpty(s)?null:s;
    		}
    	
    	void set(STRING_PROPS p,String s)
			{
			this._strings[p.ordinal()]=isEmpty(s)?".":s;
			}
    	
    	void set(INT_PROPS p,int v)
			{
			this._ints[p.ordinal()]=v<0?-1:v;
			}

    	}
   
    
	public BamStats02()
		{
		
		}
	
	
	
	private void run(String filename,SamReader r,PrintWriter out)
		{
		final Counter<Category> counter=new Counter<>();
		SAMRecordIterator iter=null;
		try
			{
			iter=r.iterator();
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(r.getFileHeader().getSequenceDictionary());
			while(iter.hasNext())
				{
				final SAMRecord record=progress.watch(iter.next());
				final Category cat=new Category();
				cat.set(STRING_PROPS.filename,filename);
				cat.flag = record.getFlags();
				cat.set(INT_PROPS.inTarget,-1);
				
				final SAMReadGroupRecord g=record.getReadGroup();
				if(g!=null) {
				cat.set(STRING_PROPS.samplename,g.getSample());
				cat.set(STRING_PROPS.platform,g.getPlatform());
				cat.set(STRING_PROPS.platformUnit,g.getPlatformUnit());
				cat.set(STRING_PROPS.library,g.getLibrary());
				}
				
				final ShortReadName readName=ShortReadName.parse(record);
				if(readName.isValid())
					{
					cat.set(STRING_PROPS.instrument,readName.getInstrumentName());
					cat.set(STRING_PROPS.flowcell,readName.getFlowCellId());
					cat.set(INT_PROPS.lane,readName.getFlowCellLane());
					cat.set(INT_PROPS.run,readName.getRunId());
					}
				
				
				if(record.getReadPairedFlag() && !record.getMateUnmappedFlag())
					{
					cat.set(STRING_PROPS.mate_chromosome,record.getMateReferenceName());
					}
				
				if(!record.getReadUnmappedFlag())
					{
					cat.set(INT_PROPS.mapq,(int)(Math.ceil(record.getMappingQuality()/10.0)*10));
					cat.set(STRING_PROPS.chromosome,record.getReferenceName());
					
					if(this.intervals!=null)
						{
						if(this.intervals.containsOverlapping(
								new Interval(record.getReferenceName(), record.getAlignmentStart(), record.getAlignmentEnd())
								))
								{

								cat.set(INT_PROPS.inTarget,1);
								}
						else
								{
								cat.set(INT_PROPS.inTarget,0);
								}
						}
					
					}
				counter.incr(cat);
				}
			progress.finish();
			
			for(final Category cat:counter.keySetDecreasing())
				{
				cat.print(out);
	    		out.print("\t");
	    		out.print(counter.count(cat));
				out.println();
				}
		
			out.flush();
			}
		finally
			{
			CloserUtil.close(iter);
			}
		
		}
	
	@Override
	public Collection<Throwable> call() throws Exception {
		final List<String> args = getInputFiles();
		SamReader samFileReader=null;
		PrintWriter out=null;
		try
			{
			if(bedFile!=null)
				{
				LOG.info("Reading BED file "+bedFile);
				this.intervals= super.readBedFileAsBooleanIntervalTreeMap(bedFile);
				}
			out = 	super.openFileOrStdoutAsPrintWriter();
			boolean first=true;
			out.print("#");
			for(final STRING_PROPS p:STRING_PROPS.values())
				{
				if(!first) out.print("\t");
				first=false;
				out.print(p.name());
				}
			for(final INT_PROPS p:INT_PROPS.values())
				{
				out.print("\t");
				out.print(p.name());
				}
    		for(final SAMFlag flg:SAMFlag.values())
    			{
    			out.print("\t");
    			out.print(flg.name());
    			}
    		out.print("\t");
    		out.print("count");
			out.println();

			
			final SamReaderFactory srf= super.createSamReaderFactory();
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				samFileReader= srf.open(SamInputResource.of(stdin()));
				run("stdin",samFileReader,out);
				samFileReader.close();
				samFileReader=null;
				}
			else
				{
				for(final String filename:IOUtils.unrollFiles(args))
					{
					LOG.info("Reading from "+filename);
					samFileReader=srf.open(new File(filename));
					run(filename,samFileReader,out);
					samFileReader.close();
					samFileReader=null;
					}
				}
			
			out.flush();
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(samFileReader);
			CloserUtil.close(out);
			}
		}
	
	public static void main(String[] args)
		{
		new BamStats02().instanceMainWithExit(args);
		}
	}
