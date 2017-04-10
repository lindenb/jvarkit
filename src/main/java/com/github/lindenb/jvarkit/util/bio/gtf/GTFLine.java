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
	public void setContig(String contig) {
		this.contig = contig;
	}
	public String getSource() {
		return source;
	}
	public void setSource(String source) {
		this.source = source;
	}
	public String getType() {
		return type;
	}
	public void setType(String type) {
		this.type = type;
	}
	public int getStart() {
		return start;
	}
	public void setStart(int start) {
		this.start = start;
	}
	public int getEnd() {
		return end;
	}
	public void setEnd(int end) {
		this.end = end;
	}
	public Double getScore() {
		return score;
	}
	public void setScore(Double score) {
		this.score = score;
	}
	public char getStrand() {
		return strand;
	}
	public void setStrand(char strand) {
		this.strand = strand;
	}
	public Map<String, String> getAtts() {
		return atts;
	}
	public void setAtts(Map<String, String> atts) {
		this.atts = atts;
	}
	public int getPhase() {
		return phase;
	}
	public void setPhase(int phase) {
		this.phase = phase;
	}
	
	private String contig;
	private String source;
	private String type;
	private int start;
	private int end;
	private Double score;
	private int phase;
	private char strand;
	private Map<String, String> atts=null;
	
}
