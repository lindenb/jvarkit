package com.github.lindenb.jvarkit.util.bio.gtf;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;

import java.util.Map;

public class GTFLine implements Locatable
	{
	public GTFLine()
		{
		
		}
	
	public Interval getInterval() {
		return new Interval(getContig(), getStart(), getEnd());
		}
	public String getContig() {
		return contig;
	}
	
	public String getSource() {
		return source;
	}
	
	public String getType() {
		return type;
	}
	
	public int getStart() {
		return start;
	}
	
	public int getEnd() {
		return end;
	}

	public Double getScore() {
		return score;
	}

	public char getStrand() {
		return strand;
	}

	public Map<String, String> getAtts() {
		return atts;
	}

	public int getPhase() {
		return phase;
	}
	
	 String contig;
	 String source;
	 String type;
	 int start;
	 int end;
	 Double score;
	 int phase;
	 char strand;
	 Map<String, String> atts=null;
	
	@Override
	public String toString() {
		final StringBuilder sb=new StringBuilder();
		sb.append(contig).append("\t").
			append(source).append("\t").
			append(type).append("\t").
			append(start).append("\t").
			append(end).append("\t").
			append(score).append("\t").
			append(phase).append("\t").
			append(strand).append("\t").
			append(atts);
		return sb.toString();
		}
}
