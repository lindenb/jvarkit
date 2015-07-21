package com.github.lindenb.jvarkit.util.bio.gtf;


import java.io.BufferedReader;
import java.io.IOException;
import java.net.URLDecoder;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;


public class GTFReader {
	private static final String GFF_VERSION="##gff-version";
	private GTFCodec codec = new GTFCodec3();
	public GTFReader() {
		}
	
	public GTFLine next(BufferedReader r) throws IOException
		{
		String line;
		while((line=r.readLine())!=null)
			{
			if(line.isEmpty()) continue;
			if(line.startsWith("#"))
				{
				if(line.startsWith(GFF_VERSION+" "))
					{
					String version =line.substring(GFF_VERSION.length()).trim();
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
				return codec.parse(line);
				}
			}
		return null;
		}
	
	
	
	private abstract class GTFCodec
		{
		Pattern tab=Pattern.compile("[\t]");
		Pattern semicolon=Pattern.compile("[;]");
		abstract GTFLine parse(String line)throws IOException;
		}
	private class GTFCodec3 extends GTFCodec
		{
		private String unescape(String s) throws IOException
			{
			return URLDecoder.decode(s, "UTF-8");
			}
		@Override
		GTFLine parse(String line) throws IOException
			{
			if(line==null) return null;
			String tokens[]=this.tab.split(line);
			GTFLine L=new GTFLine();
			L.setContig(tokens[0]);
			L.setSource(tokens[1]);
			L.setType(tokens[2]);
			L.setStart(Integer.parseInt(tokens[3]));
			L.setEnd(Integer.parseInt(tokens[4]));
			if(!tokens[5].equals(".")) L.setScore(Double.parseDouble(tokens[5]));
			L.setStrand(tokens[6].charAt(0));
			if(!tokens[7].equals(".")) L.setPhase(Integer.parseInt(tokens[7]));
			String atts[]  = semicolon.split(tokens[8]);
			Map<String, String> attMap=new HashMap<>(atts.length);
			for(String att: atts)
				{
				if(att.isEmpty()) continue;
				int eq=att.indexOf("=");
				if(eq<=0) throw new IOException("bad att "+att+" in "+line);
				String key = att.substring(0,eq);
				if(attMap.containsKey(key))  throw new IOException("duplicate "+key+" in "+line);
				String value = att.substring(eq+1);
				attMap.put(key, unescape(value) );
				}
			L.setAtts(attMap);
			return L;
			}
		}
	}
