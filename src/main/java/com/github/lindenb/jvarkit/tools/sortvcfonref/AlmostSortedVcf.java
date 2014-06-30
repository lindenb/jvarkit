package com.github.lindenb.jvarkit.tools.sortvcfonref;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import htsjdk.tribble.readers.LineIterator;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;

import htsjdk.samtools.util.CloserUtil;

/**
 * AlmostSortedVcfOnRef
 *
 */
public class AlmostSortedVcf extends AbstractCommandLineProgram
	{
	private PrintWriter out=null;
    private int MAX_RECORDS_IN_RAM=1000;
    
    @Override
    protected String getOnlineDocUrl() {
    	return "https://github.com/lindenb/jvarkit/wiki/AlmostSortedVcf";
    	}
	
    private AlmostSortedVcf()
    	{
    	
    	}
    
    @Override
    public String getProgramDescription() {
    	return "Sort an 'almost' sorted VCF. Most variants should be sorted but a few  consecutive lines might have been switched by a caller.";
    	}
    private static class ChromPosLine
    	implements Comparable<ChromPosLine>
    	{
    	int pos;
    	String chrom;
    	String line;
    	public ChromPosLine(String line)
    		{
    		int tab1=line.indexOf('\t');
    		if(tab1==-1) throw new IllegalArgumentException("Bad VCF line in "+line);
    		this.chrom=line.substring(0,tab1);
			int tab2=line.indexOf('\t',tab1+1);
    		if(tab2==-1) throw new IllegalArgumentException("Bad VCF line in "+line);
    		try
    			{
    			this.pos=Integer.parseInt(line.substring(tab1+1, tab2));
    			}
    		catch(NumberFormatException err)
    			{
    			throw new IllegalArgumentException("Bad VCF line in "+line);
    			}
    		
    		this.line=line;
    		}
    	@Override
    	public int compareTo(ChromPosLine o)
    		{
    		int i=this.chrom.compareTo(o.chrom);
    		if(i!=0) return i;
    		i=this.pos-o.pos;
    		if(i!=0) return i;
    		return this.line.compareTo(o.line);
    		}
    	@Override
    	public int hashCode() {
    		return line.hashCode();
    		}
    	@Override
    	public boolean equals(Object obj) {
    		return line.equals(ChromPosLine.class.cast(obj).line);
    		}
    	@Override
    	public String toString() {
    		return line;
    		}
    	}
	
	
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -N (int) max records in ram. default: "+MAX_RECORDS_IN_RAM);
		super.printOptions(out);
		}
    
    @Override
    public int doWork(String[] args)
    	{
    
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"N:"))!=-1)
			{
			switch(c)
				{
				case 'N': MAX_RECORDS_IN_RAM=Math.max( Integer.parseInt(opt.getOptArg()),10);break;
				default:
					{
					switch(handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
    	
		LineIterator in=null;
		int n_fixed=0;
		try
			{
			this.out=new PrintWriter(System.out);
		    ArrayList<ChromPosLine> buffer=new ArrayList<ChromPosLine>(this.MAX_RECORDS_IN_RAM);
			int ret=0;
			if(opt.getOptInd()==args.length)
				{
				in=IOUtils.openStdinForLineIterator();
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				in=IOUtils.openURIForLineIterator(filename);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			
			while(in.hasNext() && in.peek().startsWith("#"))
				{
				out.println(in.next());
				}
			
			while(in.hasNext() && !out.checkError())
				{
				String line=in.next();
				if(line.startsWith("#"))
					{
					throw new IOException("Bad VCF file in "+line);
					}
				ChromPosLine cpl=new ChromPosLine(line);
				
				if(buffer.isEmpty() )
					{
					buffer.add(cpl);
					}
				else
					{
					//not same chrom buffer : flush buffer ?
					if(!buffer.get(0).chrom.equals(cpl.chrom))
						{
						for(ChromPosLine other: buffer)
							{
							this.out.println(other.line);
							}
						buffer.clear();
						buffer.add(cpl);
						}
					else //find buffer index
						{
						int index= buffer.size();
						while(index>0 && buffer.get(index-1).pos> cpl.pos)
							{
							index--;
							}
						if(index==0) throw new IOException("Too many variants misplaced. use a larger 'max-records-in-ram' or use a classic vcf-sorter");
						if(index!=buffer.size()) ++n_fixed;
						buffer.add(index, cpl);
						
						//flush buffer if needed
						while(!buffer.isEmpty() &&buffer.size() > MAX_RECORDS_IN_RAM)
							{
							this.out.println(buffer.remove(0).line);
							}

						}
					}
				
				
				}
			
			//flush buffer
			while(!buffer.isEmpty())
				{
				this.out.println(buffer.remove(0).line);
				}
			
			this.out.flush();
			if(n_fixed==0)
				{
				info("No VCF position was fixed: worth using "+getProgramName()+" !");
				}
			else if(n_fixed>0)
				{
				info(""+n_fixed+" VCF positions were fixed: worth using "+getProgramName()+ " !");
				}
			return ret;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			CloserUtil.close(in);
			}
    	}
    
	/**
	 * @param args
	 */
	public static void main(String[] args) {
	new AlmostSortedVcf().instanceMainWithExit(args);

	}

}
