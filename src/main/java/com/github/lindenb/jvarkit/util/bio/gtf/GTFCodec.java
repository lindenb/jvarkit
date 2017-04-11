package com.github.lindenb.jvarkit.util.bio.gtf;


import java.io.IOException;
import java.net.URLDecoder;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.lang.JvarkitException;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.tribble.readers.LineIterator;


public class GTFCodec {
	private final Pattern tab=Pattern.compile("[\t]");

	private static final String GFF_VERSION="##gff-version";
	private GTFHeader header=null;

	public static class GTFHeader
		{
		private final List<String> lines = new ArrayList<>();
		private boolean is_gff3=false;
		public boolean isGff3() {
			return is_gff3;	
			}
		}

	
	
	public GTFCodec() {
		}
	
	
	
	
	public GTFLine decode(LineIterator r)throws IOException
		{
		for(;;)
			{
			if(!r.hasNext()) return null;
			String line=r.next();
			if(line.startsWith("#")) continue;
			return decode(line);	
			}
		}
		
	
	public GTFHeader readHeader(LineIterator r) throws IOException {
		if(this.header!=null) throw new IOException("Reader already read");
		this.header = new GTFHeader();
		while(r.hasNext() && r.peek().startsWith("#"))
			{
			final String line=r.next();
			if(line.startsWith(GFF_VERSION+" "))
				{
				final String version =line.substring(GFF_VERSION.length()).trim();
				if(version.equals("3"))
					{
					this.header.is_gff3=true;
					}
				}
			this.header.lines.add(line);
			}
		return this.header;
		}
		
	private String unescape(String s) throws IOException
		{
		return URLDecoder.decode(s, "UTF-8");
		}
		
	public IntervalTreeMap<List<GTFGene>> readAllAsIntervalTreeMap(final LineIterator iter) throws IOException {
		IntervalTreeMap<List<GTFGene>> h=new IntervalTreeMap<>();
		for(;;)
			{
			final List<GTFGene> L2 = nextGenesInContig(iter);
			if(L2.isEmpty()) break;
			for(final GTFGene gene:L2)
				{
				final Interval interval =new Interval(gene.getContig(),gene.getStart(),gene.getEnd());
				List<GTFGene> x= h.get(interval);
				if(x==null) {
					x=new ArrayList<>();
					h.put(interval, x);
					}
				x.add(gene);
				}
			}
		return h;
		}
		
	public List<GTFGene> readAll(final LineIterator iter) throws IOException {
		final List<GTFGene> L = new ArrayList<>();
		for(;;)
			{
			final List<GTFGene> L2 = nextGenesInContig(iter);
			if(L2.isEmpty()) break;
			L.addAll(L2);
			}
		return L;
		}
			
	public List<GTFGene> nextGenesInContig(final LineIterator iter) throws IOException {
		if(!iter.hasNext()) return Collections.emptyList();
		final Map<String,List<GTFLine>> transcript2map = new HashMap<>();
		String prevContig = null;
		while(iter.hasNext())
			{
			final String line = iter.peek();
			final GTFLine record = this.decode(line);
			if(prevContig!=null && prevContig.equals(record.getContig())) {
				break;
				}
			iter.next();//consumme
			final String transcript_id = record.getAtts().get("transcript_id");
			if( transcript_id == null ) continue;
			List<GTFLine> lines = transcript2map.get(transcript_id);
			if( lines ==null ) {
				lines = new ArrayList<>();
				transcript2map.put(transcript_id,lines);
				}
			lines.add(record);
			prevContig = record.getContig();
			}
		
		final List<GTFGene> list = new ArrayList<>(transcript2map.size());
		Collections.sort(list, (A,B)->{
			int i=A.getContig().compareTo(B.getContig());
			if( i != 0 ) return i;
			i = A.getStart() - B.getStart();
			if( i !=0 ) return i;
			i = A.getEnd() - B.getEnd();
			return i;
			});
		return list;
		}
		
	public GTFLine decode(final String line) throws IOException
		{
		if(this.header==null) {
			throw new IOException("header was not parsed");
		}
		final String tokens[]=this.tab.split(line);
		if(tokens.length<8)
			{	
			throw new JvarkitException.TokenErrors("Expected 8 columns",tokens);
			}
		final GTFLine L=new GTFLine();
		L.contig=tokens[0];
		L.source= tokens[1];
		L.type = tokens[2];
		L.start = Integer.parseInt(tokens[3]);
		L.end = Integer.parseInt(tokens[4]);
		if(!tokens[5].equals(".")) L.score = (Double.parseDouble(tokens[5]));
		L.strand = (tokens[6].charAt(0));
		if(!tokens[7].equals(".")) L.phase=Integer.parseInt(tokens[7]);
		final Map<String, String> attMap=new HashMap<>();
		final String mapStr = tokens[8];
		int k=0;
		while( k < mapStr.length())
			{
			if(Character.isWhitespace(mapStr.charAt(k))) {
				++k;
				continue;
				}
			char c= mapStr.charAt(k);
			if(c==';') { ++k; continue;}
			/* read KEY */
			final StringBuilder sbk=new StringBuilder();
			while( k < mapStr.length()) {
				c= mapStr.charAt(k);
				++k;
				if(c=='=' || Character.isWhitespace(c))
					{
					break;
					}
				sbk.append(c);
				}
			/* SKIP WS */
			while( k < mapStr.length() && Character.isWhitespace(mapStr.charAt(k))) {
				++k;
				continue;
				}
			/* EQUAL SIGN */
			if( k < mapStr.length() && mapStr.charAt(k)=='=') {
				++k;
				}
			/* SKIP WS */
			while( k < mapStr.length() && Character.isWhitespace(mapStr.charAt(k))) {
				++k;
				continue;
				}
			/* read VALUE */
			final StringBuilder sbv=new StringBuilder();
			c=(k < mapStr.length()?mapStr.charAt(k):'\0');
			// quoted string
			if( c == '\"')
				{
				++k;
				while( k < mapStr.length()) {
					c= mapStr.charAt(k);
					++k;
					if(c=='\\')
						{
						c=(k < mapStr.length()?mapStr.charAt(k):'\0');
						++k;
						switch(c) {
							case '"': sbv.append("\"");break;
							case '\'': sbv.append("\'");break;
							case 't': sbv.append("\t");break;
							case 'n': sbv.append("\n");break;
							default:break;
							}
						}
					else if(c=='\"')
						{
						break;
						}
					else
						{
						sbv.append(c);
						}
					}
				}
			else
				{
				while( k < mapStr.length()) {
					c= mapStr.charAt(k);
					++k;
					if(c==';' || Character.isWhitespace(c))
						{
						break;
						}
					sbv.append(c);
					}
				}
			attMap.put(sbk.toString(),sbv.toString());
			}
		L.atts=attMap;
		return L;
		}
	}
