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
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.Comparator;
import java.util.List;

import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

/**

BEGIN_DOC

## Example:

```bash
java -jar dist/biostar84786.jar  /path/to/input.tsv 
```

END_DOC
 */
@Program(name="biostar84786",
	biostars=84786,
	description="Matrix transposition",
	keywords={"matrix","util"}
	)
public class Biostar84786 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar84786.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names="-d",description="column delimiter")
	private char delim='\t';
	
			
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
	
	
	private int doWork(String filename,final String DELIM,PrintWriter pw) {
		if(DELIM.length()!=1)
			{
			LOG.error("DELIM must have length==1 . Got "+DELIM.length());
			return -1;
			}
		InputStream in=null;
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
			if(filename!=null)
				{
				LOG.info("opening "+filename);
				in=IOUtils.openURIForReading(filename);
				}
			else
				{
				in = stdin();
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
			if(filename!=null) in.close();
			in=null;
			CloseableIterator<Cell> iter=sorter.iterator();
			long curr_col=-1L;
			long x=0L;
			for(;;)
				{
				
				if(!iter.hasNext())
					{
					pw.println();
					break;
					}
				Cell c=iter.next();
				if(c.col!=curr_col)
					{
					if(curr_col!=-1L) pw.println();
					x=0L;
					curr_col=c.col;
					}
				if(x>0L) pw.print(DELIM);
				pw.print(c.content);
				x++;
				}
			iter.close();
			pw.flush();
			LOG.info("Done.");
			} 
		catch (Exception e)
			{
			LOG.error(e);
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
	
	@Override
	public int doWork(List<String> args) {
		String delim="\t";
		PrintWriter pw=null;
		try
			{
			pw= openFileOrStdoutAsPrintWriter(outputFile);
			int ret= doWork(oneFileOrNull(args),delim,pw);
			pw.flush();
			pw.close();pw=null;
					return ret;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(pw);	
			}
		}
	

}
