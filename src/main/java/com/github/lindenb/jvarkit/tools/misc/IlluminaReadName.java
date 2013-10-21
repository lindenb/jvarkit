package com.github.lindenb.jvarkit.tools.misc;

public class IlluminaReadName
	{
	public enum Version {V1_4};
	private String instrument;
	private int runId;
	private String flowCell;
	private int lane;
	private int tile;
	private int x;
	private int y;
	private short firstInPair=-1;
	private boolean fails_filter=false;
	private int flag;
	private String index;
	
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
	
	public static IlluminaReadName parse(Version v,String s)
		{
		switch(v)
			{
			case V1_4:  return parse1_4(s);
			default: return null;
			}
		}

	
	//@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
	public static IlluminaReadName parse1_4(String s)
		{
		int i0=(s.startsWith("@")?1:0);
		int i1=1;
		IlluminaReadName irn=new IlluminaReadName();
		i1=s.indexOf(':',i0);
		if(i1==-1) return null;
		irn.instrument=s.substring(i0, i1);
		
		
		i0=i1+1;
		i1=s.indexOf(':',i0);
		if(i1==-1) return null;
		irn.runId=Integer.parseInt(s.substring(i0, i1));
		
		i0=i1+1;
		i1=s.indexOf(':',i0);
		if(i1==-1) return null;
		irn.flowCell=s.substring(i0, i1);

		
		i0=i1+1;
		i1=s.indexOf(':',i0);
		if(i1==-1) return null;
		irn.lane=Integer.parseInt(s.substring(i0, i1));
		
		i0=i1+1;
		i1=s.indexOf(':',i0);
		if(i1==-1) return null;
		irn.tile=Integer.parseInt(s.substring(i0, i1));
		
		i0=i1+1;
		i1=s.indexOf(':',i0);
		if(i1==-1) return null;
		irn.x=Integer.parseInt(s.substring(i0, i1));
		
		i0=i1+1;
		i1=s.indexOf(' ',i0);
		if(i1==-1) return null;
		irn.y=Integer.parseInt(s.substring(i0, i1));
		
		
		i0=i1+1;
		switch(s.charAt(i0))
			{
			case '1':irn.firstInPair=1; break;
			case '2':irn.firstInPair=2; break;
			default: break;
			}
		i1=i0+1;
		
		i0=i1+1;
		switch(s.charAt(i0))
			{
			case 'Y':irn.fails_filter=true; break;
			case 'N':irn.fails_filter=false; break;
			default: break;
			}
		i1=i0+1;
		
		i0=i1+1;
		i1=s.indexOf(':',i0);
		if(i1==-1) return null;
		irn.flag=Integer.parseInt(s.substring(i0, i1));
		
		i0=i1+1;
		irn.index=s.substring(i0);
		return irn;
		}
	
	@Override
	public String toString() {
		return ""+instrument+":"+runId+":"+flowCell+":"+lane+":"+tile+":"+x+":"+y+" "
					+(firstInPair==1?"1":"2")+":"+(fails_filter?'Y':'N')+":"+flag+":"+index
					;
		}
	

	

}
