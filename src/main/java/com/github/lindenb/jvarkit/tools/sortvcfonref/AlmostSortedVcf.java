package com.github.lindenb.jvarkit.tools.sortvcfonref;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;

import htsjdk.tribble.readers.LineIterator;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;

import htsjdk.samtools.util.CloserUtil;

/**
 * AlmostSortedVcfOnRef
 *
 */
public class AlmostSortedVcf extends AbstractAlmostSortedVcf
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(AlmostSortedVcf.class);

	@Override
	public Command createCommand() {
		return new MyCommand();
		}

	static private class MyCommand extends AbstractAlmostSortedVcf.AbstractAlmostSortedVcfCommand
		{    
	
	
  
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
	protected Collection<Throwable> call(String inputName) throws Exception {
		PrintWriter out=null;
		LineIterator in=null;
		int n_fixed=0;
		try
			{
			out= openFileOrStdoutAsPrintWriter();
		    ArrayList<ChromPosLine> buffer=new ArrayList<ChromPosLine>(this.MAX_RECORDS_IN_RAM);
			if(inputName==null)
				{
				in=IOUtils.openStreamForLineIterator(stdin());
				}
			else
				{
				in=IOUtils.openURIForLineIterator(inputName);
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
					return wrapException("Bad VCF file in "+line);
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
							out.println(other.line);
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
						if(index==0) return wrapException("Too many variants misplaced. use a larger 'max-records-in-ram' or use a classic vcf-sorter");
						if(index!=buffer.size()) ++n_fixed;
						buffer.add(index, cpl);
						
						//flush buffer if needed
						while(!buffer.isEmpty() &&buffer.size() > MAX_RECORDS_IN_RAM)
							{
							out.println(buffer.remove(0).line);
							}

						}
					}
				
				
				}
			
			//flush buffer
			while(!buffer.isEmpty())
				{
				out.println(buffer.remove(0).line);
				}
			
			out.flush();
			if(n_fixed==0)
				{
				LOG.info("No VCF position was fixed: worth using "+getName()+" !");
				}
			else if(n_fixed>0)
				{
				LOG.info(""+n_fixed+" VCF positions were fixed: worth using "+getName()+ " !");
				}
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(out);
			CloserUtil.close(in);
			}
    	}
		}
    
	/**
	 * @param args
	 */
	public static void main(String[] args) {
	new AlmostSortedVcf().instanceMainWithExit(args);
	}

}
