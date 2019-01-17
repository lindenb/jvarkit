/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
* 2015: watch return SAMRecord and VariantContext

*/
package com.github.lindenb.jvarkit.util.picard;

import java.io.IOException;


import com.github.lindenb.jvarkit.util.log.Logger;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFHeader;

public class SAMSequenceDictionaryProgress
	{
    private static final Logger _LOG = Logger.build(SAMSequenceDictionaryProgress.class).make();

    private Logger log = _LOG;
	private long start_ticks=-1L;
	private long last_ticks=-1L;
	private long curr_ticks=-1L;
	
	private long count=0L;
	private long print_every_n_seconds=10L;
	private String logPrefix=null;
	private Handler handler=null;
	
	private Logger getLogger() {
		return log==null?_LOG:log;
		}
	
	public SAMSequenceDictionaryProgress logger(final Logger log) {
		this.log = log;
		return this;
	}
	
	private abstract class Handler
		{
		String lastSeen=null;

		public abstract void watch(int tid,int pos);
		public abstract void watch(String chrom,int pos);
		
		
		private String speed()
			{
			return " Speed: "+
					(count/(float)(curr_ticks-start_ticks))+" record/millisec.";
			}
		
		void printVoid()
			{
			getLogger().info( (logPrefix==null?"":"["+logPrefix+"]")+" Count:"+
					count+" Elapsed: "+duration(curr_ticks-start_ticks)+
					(lastSeen==null?"":" Last: "+lastSeen)+speed()
					);
			}
		
		}
	private class WithoutDict extends Handler
		{
		@Override
		public void watch(int tid, int pos)
			{
			super.lastSeen=(tid>=0 && pos>=0?"tid="+tid+":"+pos:null);
			printVoid();
			}
		@Override
		public void watch(String chrom, int pos)
			{
			super.lastSeen=(chrom!=null && pos>=0?chrom+":"+pos:null);
			printVoid();
			}
		}
	
	private abstract class WithDict extends Handler
		{
		SAMSequenceDictionary samSequenceDictionary;
		WithDict(SAMSequenceDictionary samSequenceDictionary)
			{
			this.samSequenceDictionary=samSequenceDictionary;
			}
		@Override
		public void watch(String chrom, int pos)
			{
			if(chrom==null || pos<0)
				{
				super.lastSeen=null;
				printVoid();
				return ;
				}
			SAMSequenceRecord ssr=this.samSequenceDictionary.getSequence(chrom);
			if(ssr==null || pos>=ssr.getSequenceLength())
				{
				super.lastSeen=null;
				printVoid();
				return ;
				}
			watch(ssr.getSequenceIndex(),pos);
			}
		
		}
	
	private  class WithOrderedDict extends WithDict
		{
		private int prev_tid=-1;
		private int prev_pos=-1;

		private long cumulLengthDone[];
		private long referenceLength=0L;
		
		WithOrderedDict(SAMSequenceDictionary dict)
			{
			super(dict);
			this.cumulLengthDone=new long[dict.size()];
			long prev_cumul=0L;
			this.referenceLength=dict.getReferenceLength();
			for(int i=0;i< dict.size();++i)
				{
				cumulLengthDone[i]=prev_cumul;
				SAMSequenceRecord ssr=dict.getSequence(i);
				prev_cumul += ssr.getSequenceLength();
				}
			}
		
		@Override
		public void watch(int tid, int pos)
			{
			if(tid==-1 || pos<0 ||  tid>=samSequenceDictionary.size())
				{
				super.lastSeen=null;
				printVoid();
				return;
				}
			
			super.lastSeen=samSequenceDictionary.getSequence(tid).getSequenceName()+":"+pos;
				
			if(prev_tid==-1 || prev_tid!=tid)
				{
				prev_tid=tid;
				prev_pos=pos;
				}
			else if(prev_tid==tid)
				{
				if(pos<prev_pos) 
					{
					getLogger().info((logPrefix==null?"":"["+logPrefix+"]")+
						"Data are not ordered on chromosome "+
						samSequenceDictionary.getSequence(tid).getSequenceName()+
						" saw "+pos+" after "+prev_pos
						);
					SAMSequenceDictionaryProgress.this.handler=new WithUnorderedDict(this.samSequenceDictionary);
					return;
					}
				prev_pos=pos;
				}
				
			
		
			
			long numBasesDone=(tid==0?0:this.cumulLengthDone[tid-1])+pos;
			long numBasesRemains=Math.max(0,referenceLength-numBasesDone);
			
			double percentDone=numBasesDone/(double)this.referenceLength;
			double millisecPerBase=(double)(curr_ticks- start_ticks)/numBasesDone;
			long timeRemain=(long)(numBasesRemains*millisecPerBase);
			
			
			getLogger().info(
					String.format("%sCount: %d Elapsed: %s(%.2f%%) Remains: %s(%.2f%%) Last: %s:%d",
					
					(logPrefix==null?"":"["+logPrefix+"]"),
					count,
					
					duration(curr_ticks-start_ticks),
					(percentDone*100.0),
					
					duration(timeRemain),
					(100-percentDone*100.0),
					
					samSequenceDictionary.getSequence(tid).getSequenceName(),
					pos
					));
			}
		}
	
	private  class WithUnorderedDict extends WithDict
		{
		WithUnorderedDict(SAMSequenceDictionary dict)
			{
			super(dict);
			}
		
		@Override
		public void watch(int tid, int pos)
			{
			if(tid==-1 || pos<0 ||  tid>=samSequenceDictionary.size())
				{
				super.lastSeen=null;
				printVoid();
				return;
				}
			
			super.lastSeen=samSequenceDictionary.getSequence(tid).getSequenceName()+":"+pos;
			printVoid();
			}
		}
	
	/**
	 * SAMSequenceDictionayProgress
	 */
	public SAMSequenceDictionaryProgress(SAMSequenceDictionary dict)
		{
		if(dict!=null)
			{
			this.handler=new WithOrderedDict(dict);
			}
		else
			{
			this.handler=new WithoutDict();
			}
		}
	
	public SAMSequenceDictionaryProgress(VCFHeader h)
		{
		this(h==null?null:h.getSequenceDictionary());
		}
	
	public SAMSequenceDictionaryProgress(final SAMFileHeader h)
		{
		this(h==null?null:h.getSequenceDictionary());
		}
	
	public void setLogPrefix(final String logPrefix)
		{
		this.logPrefix = logPrefix;
		}
	public SAMSequenceDictionaryProgress prefix(final String p)
		{
		this.setLogPrefix(p);
		return this;
		}
	
	
	public void setPrintEveryNSeconds(long print_every_n_seconds)
		{
		this.print_every_n_seconds = print_every_n_seconds;
		}
	
	public VariantContext watch(final VariantContext ctx)
		{
		if(ctx!=null)
			{
			watch(ctx.getContig(),ctx.getStart());
			}
		return ctx;
		}
	
	public SAMRecord watch(final SAMRecord rec)
		{
		if(rec!=null)
			{
			watch(rec.getReferenceIndex(),rec.getAlignmentStart());
			}
		return rec;
		}
	
	private boolean incr()
		{
		this.count++;
		this.curr_ticks=System.currentTimeMillis();
	
		if(this.start_ticks==-1L )
			{
			this.start_ticks=curr_ticks;
			this.last_ticks=curr_ticks;
			return false;
			}
		if(curr_ticks-this.last_ticks  < print_every_n_seconds*1000) return false;
		return true;
		}
	
	public void watch(String chrom,int pos)
		{
		if(!incr()) return;
		this.handler.watch(chrom, pos);
		last_ticks=curr_ticks;
		}
	
	public void watch(int tid,int pos)
		{
		if(!incr()) return;
		this.handler.watch(tid,pos);
		last_ticks=curr_ticks;
		}	
	
	private String duration(long millisecs)
		{
		long n =millisecs/1000;
		if(n<60)
			{
			return n+" second"+(n<2?"":"s");
			}
		n/=60;//minutes
		if(n< 60)
			{
			return n+" minute"+(n<2?"":"s");
			}
		n/=60;//hours
		if(n< 24)
			{
			return n+" hour"+(n<2?"":"s");
			}
		n/=24;
	
		if(n< 365)
			{
			return n+" day"+(n<2?"":"s");
			}
		n/=365;
		
		return n+" year"+(n<2?"":"s");
		}
	
	/** return the number of records seen so far */
	public long getCount()
		{
		return this.count;
		}
	
	public void finish()
		{
		getLogger().info("done: N="+getCount());
		}

	public static VCFIterator wrap(Logger log,VCFIterator r)
		{
		return new LogVCFIterator(log,r);
		}
	
	private static class LogVCFIterator  implements VCFIterator
		{
		private final VCFIterator delegate;
		private final SAMSequenceDictionaryProgress progress;
		LogVCFIterator(final Logger log,VCFIterator delegate) {
			this.delegate=delegate;
			this.progress=new SAMSequenceDictionaryProgress(delegate.getHeader()).logger(log);
			}
		//@Override
		//public AbstractVCFCodec getCodec() { return this.delegate.getCodec(); }
		@Override
		public VCFHeader getHeader() { return this.delegate.getHeader();}
		@Override
		public VariantContext peek() { return this.delegate.peek();}

		@Override
		public void close()  {
			this.delegate.close();
			this.progress.finish();
			}
		@Override
		public boolean hasNext() {
			return this.delegate.hasNext();
		}
		@Override
		public VariantContext next() {
			return this.progress.watch(this.delegate.next());
			}
		}
	
	public static CloseableIterator<SAMRecord> wrap(Logger log,SAMFileHeader header,CloseableIterator<SAMRecord> r)
		{
		return new LogSamRecordIterator(log,header,r);
		}

private static class LogSamRecordIterator  implements CloseableIterator<SAMRecord>
	{
	private final CloseableIterator<SAMRecord> delegate;
	private final SAMSequenceDictionaryProgress progress;
	LogSamRecordIterator(final Logger log,SAMFileHeader header,CloseableIterator<SAMRecord> delegate) {
		this.delegate=delegate;
		this.progress=new SAMSequenceDictionaryProgress(header).logger(log);
		}

	@Override
	public void close(){
		this.delegate.close();
		this.progress.finish();
		}
	@Override
	public boolean hasNext() {
		return this.delegate.hasNext();
	}
	@Override
	public SAMRecord next() {
		return this.progress.watch(this.delegate.next());
		}
	}

	
	}
