package com.github.lindenb.jvarkit.tools.genome2svg.beans;

import java.nio.file.Path;
import java.nio.file.Paths;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;

public class PathBean extends AbstractBean {
	private String filename;
	
	public PathBean() {
		this("");
		}
	public PathBean(final String fn) {
		this.filename = fn;
		}
	
	public String getPath() {
		return this.filename;
		}
	public void setPath(final String f) {
		this.filename = f;
		}
	public Path asPath() {
		return Paths.get(StringUtils.assertNotBlank(getPath(),"undefined path"));
		}
	
	public SAMSequenceDictionary getSAMSequenceDictionary() {
		return null;
		}
	
	public String resolveContig(final String ctg) {
		SAMSequenceDictionary dict = getSAMSequenceDictionary();
		if(dict==null) return ctg;
		return ctg;
		}
	public Locatable resolveContig(final Locatable loc) {
		if(loc==null) return null;
		final String s=resolveContig(loc.getContig());
		if(StringUtils.isBlank(s)) return null;
		if(s.equals(loc.getContig())) return loc;
		return new SimpleInterval(s,loc.getStart(),loc.getEnd());
		}


	}
