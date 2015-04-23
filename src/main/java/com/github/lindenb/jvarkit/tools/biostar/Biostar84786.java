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
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Comparator;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

public class Biostar84786 extends AbstractCommandLineProgram
	{
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

	private Biostar84786()
		{
		
		}
	
	@Override
	protected String getOnlineDocUrl()
		{
		return DEFAULT_WIKI_PREFIX+"Biostar86363";
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -d (delim) column delimter: default = tab");
		super.printOptions(out);
		}
	
	@Override
	public String getProgramDescription() {
		return " Matrix transposition ( see  http://www.biostars.org/p/84786/ ).";
		}

	
	private int doWork(File IN,String DELIM) {
		if(DELIM.length()!=1)
			{
			error("DELIM must have length==1 . Got "+DELIM.length());
			return -1;
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
		try {
			final char delimiter=DELIM.charAt(0);
			sorter=SortingCollection.newInstance(
				Cell.class, new CellCodec(), comparator,10000
				);
			sorter.setDestructiveIteration(true);
			if(IN!=null)
				{
				this.info("opening "+IN);
				in=IOUtils.openFileForReading(IN);
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
					if(row%10000==0) this.info("row:"+row);
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
			if(IN!=null) in.close();
			in=null;
			CloseableIterator<Cell> iter=sorter.iterator();
			long curr_col=-1L;
			long x=0L;
			for(;;)
				{
				
				if(!iter.hasNext())
					{
					System.out.println();
					break;
					}
				Cell c=iter.next();
				if(c.col!=curr_col)
					{
					if(curr_col!=-1L) System.out.println();
					x=0L;
					curr_col=c.col;
					}
				if(x>0L) System.out.print(DELIM);
				System.out.print(c.content);
				x++;
				}
			iter.close();
			this.info("Done.");
			} 
		catch (Exception e)
			{
			e.printStackTrace();
			this.error(e,"BOUM");
			return -1;
			}
		finally
			{
			if(sorter!=null) sorter.cleanup();
			if(in!=null) CloserUtil.close(in);
			}
		return 0;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Biostar84786().instanceMainWithExit(args);

	}

	
	public int doWork(String[] args)
		{
		String delim="\t";
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"d:"))!=-1)
			{
			switch(c)
				{
				case 'd': delim=opt.getOptArg();break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		try
			{
			if(opt.getOptInd()==args.length)
				{
				return doWork(null,delim);
				}
			else if(opt.getOptInd()+1==args.length)
				{
				return doWork(new File(args[opt.getOptInd()]),delim);
				}
			else
				{
				error("Illegal Number of arguments");
				return -1;
				}
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		}
	

}
