package com.github.lindenb.jvarkit.util.picard;

import java.io.IOException;
import java.util.regex.Pattern;

import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
import net.sf.samtools.SAMSequenceDictionary;

import org.broad.tribble.readers.LineReader;
import com.github.lindenb.jvarkit.lang.Function;
import com.github.lindenb.jvarkit.util.biomart.BiomartQuery;

public class IntervalTreeMapFactory<T>
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
	
	
	
	public IntervalTreeMap<T> createIntervalMap(LineReader r) throws IOException
	 	{
	 	Pattern tab=Pattern.compile("[\t]");
	 	String line;
	 	IntervalTreeMap<T> map=new IntervalTreeMap<T>();
	 	while((line=r.readLine())!=null)
             {
             String tokens[]=tab.split(line);
             Interval interval=getKeyFunction().apply(tokens);
             if(interval==null) continue;
             if(getSamSequenceDictionary()!=null &&
            		getSamSequenceDictionary().getSequence(interval.getSequence())==null 
            		) continue;
             T value=getValueFunction().apply(tokens);
             if(value==null) continue;
             map.put(interval,value);
             }
	 	return map;
	 	}
		 
	 public static void main(String[] args) throws Exception
		{
		BiomartQuery q=new BiomartQuery();
		q.setAttributes("chromosome_name","start_position", "end_position",
                "entrezgene");
		LineReader r= q.execute();
		String line;
		System.err.println("LOading");
	     int n=0;
	     IntervalTreeMap<String> map=new IntervalTreeMap<String>();
	     while((line=r.readLine())!=null && (n<100))
	             {
	             String tokens[]=line.split("\t");
	             if(tokens.length<4 || tokens[3].isEmpty())
	                     {
	                     continue;
	                     }
	             ++n;
	             map.put(new Interval(tokens[0], Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2])), line);
	             }
	     r.close();
	     System.err.println("Done." + n);
	
	     System.err.println(map.getOverlapping(new Interval("15",25453300,25453305)));//100033603        15      25453230        25453310

		}
	}
