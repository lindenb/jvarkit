package com.github.lindenb.jvarkit.util.picard;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.SAMSequenceDictionary;

import com.github.lindenb.jvarkit.lang.Function;



public abstract class AbstractIntervalMapFactory<T>
	{
	
	private Function<String[],T> funValue=null;
	private Function<String[],Interval> funKey=new DefaultIntervalFunction();
	private SAMSequenceDictionary samSequenceDictionary=null;
	
	public static class DefaultIntervalFunction
		implements Function<String[],Interval>
		{
		private int chromCol=0;
		private int startCol=1;
		private int endCol=2;
		public DefaultIntervalFunction()
			{
			this(0,1,2);
			}
		
		public DefaultIntervalFunction(int chromCol,int startCol,int endCol)
			{
			this.chromCol=chromCol;
			this.startCol=startCol;
			this.endCol=endCol;
			}
		@Override
		public Interval apply(String[] param)
			{
			if(	chromCol>=param.length ||
				startCol>=param.length ||
				endCol>=param.length) return null;
			
			return new Interval(
					param[0],
					Integer.parseInt(param[1]),
					Integer.parseInt(param[2])
					);
			}
		}
	
	protected AbstractIntervalMapFactory()
		{
		
		}

	
	public void setValueFunction(Function<String[],T> fun)
		{
		this.funValue=fun;
		}
	
	public Function<String[],T> getValueFunction()
		{
		return this.funValue;
		}
	
	public void setKeyFunction(Function<String[],Interval> fun)
		{
		this.funKey=fun;
		}
	
	public Function<String[],Interval> getKeyFunction()
		{
		return this.funKey;
		}
	
	public SAMSequenceDictionary getSamSequenceDictionary()
		{
		return samSequenceDictionary;
		}
	
	public void setSamSequenceDictionary(
			SAMSequenceDictionary samSequenceDictionary)
		{
		this.samSequenceDictionary = samSequenceDictionary;
		}
	
	protected int compareIntervals(Interval i1,Interval i2)
		{
		int i=compareChromosomes(i1.getContig(),i2.getContig());
		if(i!=0) return i;
		i=i1.getStart()-i2.getStart();
		if(i!=0) return i;
		i=i1.getEnd()-i2.getEnd();
		return i;
		}

	
	protected int compareChromosomes(String c1,String c2)
		{
		SAMSequenceDictionary d;
		if((d=getSamSequenceDictionary())!=null)
			{
			return d.getSequenceIndex(c1)-d.getSequenceIndex(c2);
			}
		return c1.compareTo(c2);
		}
	
	
	}
