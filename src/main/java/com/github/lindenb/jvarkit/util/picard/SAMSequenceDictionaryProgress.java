package com.github.lindenb.jvarkit.util.picard;

import java.util.logging.Logger;


import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

public class SAMSequenceDictionaryProgress
	{
	private static final Logger LOG=Logger.getLogger("jvarkit");

	private long start_ticks=-1L;
	private long last_ticks=-1L;
	private long curr_ticks=-1L;
	
	private long count=0L;
	private long print_every_n_seconds=10L;
	private String prefix=null;
	private Handler handler=null;
	
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
			LOG.info( (prefix==null?"":"["+prefix+"]")+" Count:"+
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
					LOG.info((prefix==null?"":"["+prefix+"]")+
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
			
			
			LOG.info(
					String.format("%sCount: %d Elapsed: %s(%.2f%%) Remains: %s(%.2f%%) Last: %s:%d",
					
					(prefix==null?"":"["+prefix+"]"),
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
	
	public void setLogPrefix(String prefix)
		{
		this.prefix = prefix;
		}
	
	public void setPrintEveryNSeconds(long print_every_n_seconds)
		{
		this.print_every_n_seconds = print_every_n_seconds;
		}
	
	
	public void watch(final SAMRecord rec)
		{
		if(rec==null) return;
		watch(rec.getReferenceIndex(),rec.getAlignmentStart());
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
		LOG.info("done: N="+getCount());
		}


	}
