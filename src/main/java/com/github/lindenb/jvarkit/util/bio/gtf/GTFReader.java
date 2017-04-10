package com.github.lindenb.jvarkit.util.bio.gtf;


import java.io.BufferedReader;
import java.io.IOException;
import java.net.URLDecoder;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.regex.Pattern;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.LineIterator;


public class GTFReader {
	private static final String GFF_VERSION="##gff-version";
	private GTFCodec codec = new GTFCodec3();
	public GTFReader() {
		}
	
	public GTFLine next(final BufferedReader r) throws IOException
		{
		String line;
		while((line=r.readLine())!=null)
			{
			if(line.isEmpty()) continue;
			if(line.startsWith("#"))
				{
				if(line.startsWith(GFF_VERSION+" "))
					{
					final String version =line.substring(GFF_VERSION.length()).trim();
					if(version.equals("3"))
						{
						this.codec = new GTFCodec3();
						}
					else
						{
						throw new IOException("Unsupported codec : "+line);
						}
					}
				}
			else
				{	
				return codec.decode(line);
				}
			}
		return null;
		}
	
	
	public static class GTFHeader
		{
		private final List<String> lines = new ArrayList<>();
		private boolean is_gff3=false;
		public boolean isGff3() {
			return is_gff3;	
			}
		}
	
	
	private abstract class GTFCodec
		{
		Pattern tab=Pattern.compile("[\t]");
		Pattern semicolon=Pattern.compile("[;]");
		public abstract GTFHeader readHeader(LineIterator r)throws IOException;
		public GTFLine decode(LineIterator r)throws IOException
			{
			return r.hasNext()?decode(r.next()):null;
			}
		public abstract GTFLine decode(String s)throws IOException;
		}
	
	private class GTFCodec3 extends GTFCodec
		{
		private GTFHeader header=null;
		private Pattern key_pattern=Pattern.compile("[a-zA-Z][_a-zA-Z_0-9]*");
		private Pattern data_pattern=Pattern.compile("\"(\\\\.|[^\\\\\"])*\"");
		@Override
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
		
		@Override
		public GTFLine decode(final String line) throws IOException
			{
			if(this.header==null) {
				throw new IOException("header was not parsed");
			}
			final String tokens[]=this.tab.split(line);
			final GTFLine L=new GTFLine();
			L.setContig(tokens[0]);
			L.setSource(tokens[1]);
			L.setType(tokens[2]);
			L.setStart(Integer.parseInt(tokens[3]));
			L.setEnd(Integer.parseInt(tokens[4]));
			if(!tokens[5].equals(".")) L.setScore(Double.parseDouble(tokens[5]));
			L.setStrand(tokens[6].charAt(0));
			if(!tokens[7].equals(".")) L.setPhase(Integer.parseInt(tokens[7]));
			final Map<String, String> attMap;
			if( this.header.isGff3())
				{
				attMap=new HashMap<>();
				final Scanner scanner= new Scanner(tokens[8]);
				for(;;)
					{
					final String key=scanner.findInLine(this.key_pattern);
					if(key==null) break;
					final String data=scanner.findInLine(this.data_pattern);
					if(data==null) throw new IOException("");
					attMap.put(key,unescape(data));
					String semicolon=scanner.findInLine("");
					if(semicolon==null) break;
					if(!semicolon.equals(";"))
						throw new RuntimeIOException("Boum");
					}
				scanner.close();
				}
			else
				{
				String atts[]  = this.semicolon.split(tokens[8]);
				attMap=new HashMap<>(atts.length);
				for(final String att: atts)
					{
					if(att.isEmpty()) continue;
					int eq=att.indexOf("=");
					
					if(eq<=0)
						{
						String key = att.substring(0,eq);
						if(attMap.containsKey(key))  throw new IOException("duplicate "+key+" in "+line);
						String value = att.substring(eq+1);
						attMap.put(key, unescape(value) );
						}
					}
				}
			L.setAtts(attMap);
			return L;
			}
		}
	}
