package com.github.lindenb.jvarkit.tools.biostar;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Comparator;

import com.github.lindenb.jvarkit.util.picard.cmdline.Option;
import com.github.lindenb.jvarkit.util.picard.cmdline.StandardOptionDefinitions;
import com.github.lindenb.jvarkit.util.picard.cmdline.Usage;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

public class Biostar84786 extends AbstractCommandLineProgram {
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" Matrix transposition ( see  http://www.biostars.org/p/84786/ ).";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="input stream default: stdin. ",
    		optional=true)
	public File IN=null;
    
    @Option(shortName= "D", doc="delimiter. Default is tabulation. ", optional=true)
    public String DELIM="\t";
    
	private Log LOG=Log.getInstance(Biostar84786.class);
	
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
	protected int doWork() {
		if(DELIM.length()!=1)
			{
			LOG.error("DELIM must have length==1 . Got "+DELIM.length());
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
			sorter=SortingCollection.newInstance(Cell.class, new CellCodec(), comparator, super.MAX_RECORDS_IN_RAM);
			sorter.setDestructiveIteration(true);
			if(IN!=null)
				{
				LOG.info("opening "+IN);
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
			LOG.info("Done.");
			} 
		catch (Exception e)
			{
			e.printStackTrace();
			LOG.error(e,"BOUM");
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

}
