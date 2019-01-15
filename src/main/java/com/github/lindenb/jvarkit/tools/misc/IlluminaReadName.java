/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


*/
package com.github.lindenb.jvarkit.tools.misc;

import java.util.Optional;
import java.util.function.Function;
import java.util.regex.Pattern;

public class IlluminaReadName
	{
	private static final Pattern COLON  = Pattern.compile("[\\:]");
	public enum Version {V1_4};
	private String instrument=null;
	private int runId=-1;
	private String flowCell=null;
	private int lane=-1;
	private int tile=-1;
	private int x=-1;
	private int y=-1;
	private short firstInPair=-1;
	private boolean fails_filter=false;
	private int flag=-1;
	private String index=null;
	
	
	public static class Parser implements Function<String, Optional<IlluminaReadName>>
		{
		private Version version = Version.V1_4;
		private boolean isPositiveInteger(final String s) {
			try {
			return Integer.parseInt(s)>=0;	
			}
		catch(NumberFormatException err) {
			return false;
			}
		} 
		
		private Optional<IlluminaReadName> parse_any(final String s) {
			int blank=0;
			for(blank=0;blank< s.length();blank++)
				{
				if(Character.isWhitespace(s.charAt(blank))) break;
				}
			String left=s.substring(0,blank);
			String right=blank>=s.length()?"":s.substring(blank).trim();
			
			
			String tokens[] = COLON.split(left);
			if(tokens[0].startsWith("@")) {
				tokens[0]=tokens[0].substring(1);
			}
			
			final IlluminaReadName irn=new IlluminaReadName();
			switch(tokens.length)
				{
				case 7:
					{
					if(!isPositiveInteger(tokens[1])) return Optional.empty();
					if(!isPositiveInteger(tokens[3])) return Optional.empty();
					if(!isPositiveInteger(tokens[4])) return Optional.empty();
					if(!isPositiveInteger(tokens[5])) return Optional.empty();
					if(!isPositiveInteger(tokens[6])) return Optional.empty();
					irn.instrument=tokens[0];
					irn.runId = Integer.parseInt(tokens[1]);
					irn.flowCell=tokens[2];
					irn.lane=Integer.parseInt(tokens[3]);
					irn.tile=Integer.parseInt(tokens[4]);
					irn.x=Integer.parseInt(tokens[5]);
					irn.y=Integer.parseInt(tokens[6]);
					break;
					}
				case 5:/* @HWUSI-EAS100R:6:73:941:1973#0/1 */
					{
					int hash=tokens[4].indexOf("#");
					if(hash!=-1)
						{
						String tokens4right=tokens[4].substring(hash+1);
						tokens[4]=tokens[4].substring(0, hash);
						if(tokens4right.endsWith("/1")) {
							irn.firstInPair=1;
							}
						else if(tokens4right.endsWith("/2")) {
							irn.firstInPair=2;
							}
						}
					if(!isPositiveInteger(tokens[1])) return Optional.empty();
					if(!isPositiveInteger(tokens[2])) return Optional.empty();
					if(!isPositiveInteger(tokens[3])) return Optional.empty();
					if(!isPositiveInteger(tokens[4])) return Optional.empty();
					irn.instrument=tokens[0];
					irn.lane = Integer.parseInt(tokens[1]);
					irn.tile=Integer.parseInt(tokens[2]);
					irn.x=Integer.parseInt(tokens[3]);
					irn.y=Integer.parseInt(tokens[4]);
					break;
					}
				default: return Optional.empty();
				}
			
			if(right.isEmpty()) return Optional.of(irn);
			tokens = COLON.split(right);
			
			if(tokens.length>0 && tokens[0].length()==1) {
				switch(tokens[0].charAt(0))
					{
					case '1':irn.firstInPair=1; break;
					case '2':irn.firstInPair=2; break;
					default: break;
					}
				}
			
			if(tokens.length>1 && tokens[1].length()==1) {
				switch(tokens[1].charAt(0))
					{
					case 'Y':irn.fails_filter=true; break;
					case 'N':irn.fails_filter=false; break;
					default: break;
					}
				}
			
			if(tokens.length>2 && isPositiveInteger(tokens[2])) {
				irn.flag=Integer.parseInt(tokens[2]);
			}
			
			if(tokens.length>3 ) {
				irn.index=tokens[3];
			}
			
			return Optional.of(irn);
		}
		
		public void setVersion(Version version) {
			this.version = version;
		}
		
		public Version getVersion() {
			return version;
		}
		
		@Override
		public Optional<IlluminaReadName> apply(final String s) {
			switch(getVersion())
				{
				case V1_4: return parse_any(s);
				default: throw new IllegalArgumentException("Cannot parse read name using version "+getVersion());
				}
			}
		}
	
	private IlluminaReadName()
		{
		
		}
	
	public String getInstrument()
		{
		return instrument;
		}
	
	public String getIndex()
		{
		return index;
		}
	
	public String getFlowCell()
		{
		return flowCell;
		}
	
	public int getLane() {
		return lane;
	}
	
	public int getRunId() {
		return runId;
	}
	
	public int getTile() {
		return tile;
	}
	
	public int getX() {
		return x;
	}
	
	public int getY() {
		return y;
	}
	
	@Deprecated //use Parser
	public static IlluminaReadName parse(Version v,String s)
		{
		switch(v)
			{
			case V1_4:  return parse1_4(s);
			default: return null;
			}
		}

	@Deprecated //use Parser
	//@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
	public static IlluminaReadName parse1_4(String s)
		{
		return new Parser().apply(s).orElse(null);
		}
	
	@Override
	public String toString() {
		return ""+instrument+":"+runId+":"+flowCell+":"+lane+":"+tile+":"+x+":"+y+" "
					+(firstInPair==1?"1":"2")+":"+(fails_filter?'Y':'N')+":"+flag+":"+index
					;
		}
	

	

}
