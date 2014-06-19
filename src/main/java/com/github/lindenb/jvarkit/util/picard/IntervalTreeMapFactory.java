package com.github.lindenb.jvarkit.util.picard;

import java.io.IOException;
import java.util.regex.Pattern;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.tribble.readers.LineReader;


public class IntervalTreeMapFactory<T>
	extends AbstractIntervalMapFactory<T>
	{
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
	}
