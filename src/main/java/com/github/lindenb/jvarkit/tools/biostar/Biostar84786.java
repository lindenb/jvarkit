/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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


History:
* 2015: moved to AbstractCommandLineProgram

*/
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

public class Biostar84786 extends AbstractBiostar84786
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Biostar84786.class);

	
	private static class Cell
		{
		long row;
		long col;
		String content;
		Cell()
			{
			
			}
		Cell(long row,long col,StringBuilder b)
			{
			this.row=row;
			this.col=col;
			this.content=b.toString();
			}
		@Override
		public String toString() {
			return "("+row+","+col+":"+content+")";
			}
		}
	
	private static class CellCodec
		extends AbstractDataCodec<Cell>
		{
		@Override
		public Cell decode(DataInputStream dis) throws IOException
			{
			Cell c=new Cell();
			try {
				c.row=dis.readLong();
			} catch (Exception e) {//EOF reached
				throw new IOException(e);
				}
			c.col=dis.readLong();
			c.content=dis.readUTF();
			return c;
			}
		@Override
		public void encode(DataOutputStream dos, Cell c)
				throws IOException {
			dos.writeLong(c.row);
			dos.writeLong(c.col);
			dos.writeUTF(c.content);
			}
		@Override
		public AbstractDataCodec<Cell> clone() {
			return new CellCodec() ;
			}
		}

	@Override
	public Command createCommand() {
		return new MyCommand();
		}
	
	static private class MyCommand extends AbstractBiostar84786.AbstractBiostar84786Command
		{
		@Override
		protected Collection<Throwable> call(String inputName) throws Exception
			{	
			if(super.delim.length()!=1)
				{
				return wrapException("DELIM must have length==1 . Got "+delim.length());
				}
			InputStream in=System.in;
			SortingCollection<Cell> sorter=null;
			final  Comparator<Cell> comparator=new Comparator<Biostar84786.Cell>() {
				@Override
				public int compare(final Cell o1, final Cell o2)
					{
					int i;
					i=(o1.col<o2.col?-1:o1.col>o2.col?1:0);
					if(i!=0) return i;
					i=(o1.row<o2.row?-1:o1.row>o2.row?1:0);
					if(i!=0) return i;
					return o1.content.compareTo(o2.content);
					}
				};
			PrintStream out=null;
			try {
				out= super.openFileOrStdoutAsPrintStream();
				
				final char delimiter=delim.charAt(0);
				sorter=SortingCollection.newInstance(
					Cell.class, new CellCodec(), comparator,10000
					);
				sorter.setDestructiveIteration(true);
				if(inputName==null)
					{
					in = stdin();
					}
				else
					{
					in = IOUtils.openURIForReading(inputName);
					}
				long row=0L;
				long col=0L;
				StringBuilder b=new StringBuilder();
				for(;;)
					{
					int c=in.read();
					if(c=='\n' || c==-1)
						{
						sorter.add(new Cell(row,col,b));
						row++;
						col=0;
						b.setLength(0);
						if(c==-1) break;
						if(row%10000==0) LOG.info("row:"+row);
						}
					else if(c==delimiter)
						{
						sorter.add(new Cell(row,col,b));
						b.setLength(0);
						col++;
						}
					else
						{
						b.append((char)c);
						}
					}
				sorter.doneAdding();
				CloserUtil.close(in);
				in=null;
				CloseableIterator<Cell> iter=sorter.iterator();
				long curr_col=-1L;
				long x=0L;
				for(;;)
					{
					
					if(!iter.hasNext())
						{
						out.println();
						break;
						}
					Cell c=iter.next();
					if(c.col!=curr_col)
						{
						if(curr_col!=-1L) out.println();
						x=0L;
						curr_col=c.col;
						}
					if(x>0L) out.print(super.delim);
					out.print(c.content);
					x++;
					}
				iter.close();
				LOG.info("Done.");
				return Collections.emptyList();
				} 
			catch (Exception e)
				{
				return wrapException(e);
				}
			finally
				{
				if(sorter!=null) sorter.cleanup();
				CloserUtil.close(in);
				CloserUtil.close(out);
				}
			}
		}
	
	public static void main(String[] args) {
		new Biostar84786().instanceMainWithExit(args);

	}

}
