/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.sortvcfonref;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;

import htsjdk.tribble.readers.LineIterator;

import com.github.lindenb.jvarkit.io.IOUtils;

import htsjdk.samtools.util.CloserUtil;

/**
 * AlmostSortedVcfOnRef
 *
 */
public class AlmostSortedVcf extends AbstractAlmostSortedVcf
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(AlmostSortedVcf.class);

    
    private static class ChromPosLine
    	implements Comparable<ChromPosLine>
    	{
    	final int pos;
    	final String chrom;
    	final String line;
    	public ChromPosLine(final String line)
    		{
    		final int tab1=line.indexOf('\t');
    		if(tab1==-1) throw new IllegalArgumentException("Bad VCF line in "+line);
    		this.chrom=line.substring(0,tab1);
    		final int tab2=line.indexOf('\t',tab1+1);
    		if(tab2==-1) throw new IllegalArgumentException("Bad VCF line in "+line);
    		try
    			{
    			this.pos=Integer.parseInt(line.substring(tab1+1, tab2));
    			}
    		catch(final NumberFormatException err)
    			{
    			throw new IllegalArgumentException("Bad VCF line in "+line);
    			}
    		
    		this.line=line;
    		}
    	@Override
    	public int compareTo(final  ChromPosLine o)
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
    	public boolean equals(final  Object obj) {
    		return line.equals(ChromPosLine.class.cast(obj).line);
    		}
    	@Override
    	public String toString() {
    		return line;
    		}
    	}
	
	@Override
	protected Collection<Throwable> call(final  String inputName) throws Exception {
			
		PrintWriter out=null;
		LineIterator in=null;
		int n_fixed=0;
		try
			{
			out= super.openFileOrStdoutAsPrintWriter();
			final ArrayList<ChromPosLine> buffer=new ArrayList<ChromPosLine>(this.MAX_RECORDS_IN_RAM);
			if(inputName==null)
				{
				in=IOUtils.openStdinForLineIterator();
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
				final  String line=in.next();
				if(line.startsWith("#"))
					{
					out.close();
					throw new IOException("Bad VCF file in "+line);
					}
				final  ChromPosLine cpl=new ChromPosLine(line);
				
				if(buffer.isEmpty() )
					{
					buffer.add(cpl);
					}
				else
					{
					//not same chrom buffer : flush buffer ?
					if(!buffer.get(0).chrom.equals(cpl.chrom))
						{
						for(final  ChromPosLine other: buffer)
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
						if(index==0) {
							out.close();
							throw new IOException("Too many variants misplaced. use a larger 'max-records-in-ram' or use a classic vcf-sorter");
						}
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
		catch(final  Exception err)
			{
			return wrapException(err);
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
