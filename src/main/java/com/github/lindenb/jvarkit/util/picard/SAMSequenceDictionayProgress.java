package com.github.lindenb.jvarkit.util.picard;

import java.util.logging.Logger;


import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

public class SAMSequenceDictionayProgress
	{
	private static final Logger LOG=Logger.getLogger("jvarkit");
	private SAMSequenceDictionary samSequenceDictionary;
	private long start=-1L;
	private long last=-1L;
	private long cumulLengthDone[];
	private long referenceLength=0L;
	private long count=0L;
	private long print_every_n_seconds=10L;
	private String prefix=null;
	
	/**
	 * SAMSequenceDictionayProgress
	 */
	public SAMSequenceDictionayProgress(SAMSequenceDictionary dict)
		{
		this.samSequenceDictionary=dict;
		if(dict!=null)
			{
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
		}
	
	public void setLogPrefix(String prefix)
		{
		this.prefix = prefix;
		}
	
	public void setPrintEveryNSeconds(long print_every_n_seconds)
		{
		this.print_every_n_seconds = print_every_n_seconds;
		}
	
	public SAMSequenceDictionary getSamSequenceDictionary()
		{
		return samSequenceDictionary;
		}
	
	public void watch(final SAMRecord rec)
		{
		if(rec==null) return;
		watch(rec.getReferenceIndex(),rec.getAlignmentStart());
		}
	
	
	public void watch(String chrom,int pos)
		{
		if(chrom==null || samSequenceDictionary==null) return;
		watch(samSequenceDictionary.getSequenceIndex(chrom),pos);
		}
	
	private String duration(long millisecs)
		{
		long n =millisecs/1000;
		if(n<60)
			{
			return n+" seconds";
			}
		n/=60;//minutes
		if(n< 60)
			{
			return n+" minutes";
			}
		n/=60;//hours
		if(n< 24)
			{
			return n+" hours";
			}
		n/=24;
	
		return n+" days";
		}
	
	/** return the number of records seen so far */
	public long getCount()
		{
		return this.count;
		}
	
	public void watch(int tid,int pos)
		{
		this.count++;
		long curr=System.currentTimeMillis();
	
		if(start==-1L )
			{
			start=curr;
			last=curr;
			return;
			}
		
		if(curr-last < print_every_n_seconds*1000) return;
		
		if(samSequenceDictionary==null)
			{
			
			LOG.info(
					String.format("%sCount: %,10d",
					(this.prefix==null?"":"["+this.prefix+"]"))+
					count
					);
			return;
			}

		
		if(tid<0 || tid>=samSequenceDictionary.size() || pos<1) return;
	
		
		long numBasesDone=(tid==0?0:this.cumulLengthDone[tid-1])+pos;
		long numBasesRemains=Math.max(0,referenceLength-numBasesDone);
		
		double percentDone=numBasesDone/(double)this.referenceLength;
		double millisecPerBase=(double)(curr-this.start)/numBasesDone;
		long timeRemain=(long)(numBasesRemains*millisecPerBase);
		
		
		LOG.info(
				String.format("%sCount: %d Elapsed: %s(%.2f%%) Remains: %s(%.2f%%) Last: %s:%d",
				
				(this.prefix==null?"":"["+this.prefix+"]"),
				count,
				
				duration(curr-start),
				(percentDone*100.0),
				
				duration(timeRemain),
				(100-percentDone*100.0),
				
				samSequenceDictionary.getSequence(tid).getSequenceName(),
				pos
				));
		this.last=curr;
		}

	}
