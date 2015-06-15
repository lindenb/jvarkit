package com.github.lindenb.jvarkit.util.bio.bed;

import htsjdk.tribble.Feature;

public class BedLine
	implements Feature
	{
	private String tokens[];
	
	public BedLine(String tokens[])
		{
		this.tokens=tokens;
		}
	@Override
	@Deprecated
	public final String getChr() {
		return getContig();
		}
	
	@Override
	public String getContig() {
		return tokens[0];
		}
	@Override
	public int getStart() {
		return Integer.parseInt(tokens[1]);
		}
	
	@Override
	public int getEnd() {
		return (tokens.length<3 ?getStart(): Integer.parseInt(tokens[2]));
		}

	public String get(int index)
		{
		return (index<tokens.length?tokens[index]:null);
		}
	
	public int getColumnCount()
		{
		return tokens.length;
		}
	public static boolean isBedHeader(String line)
		{
		return line.startsWith("#") || line.startsWith("track") || line.startsWith("browser");
		}

	}
