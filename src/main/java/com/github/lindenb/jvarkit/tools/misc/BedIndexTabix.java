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
* 2015 : moving to knime
* 2014 : creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Comparator;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.util.TabixUtils;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;


public class BedIndexTabix
	extends AbstractBedIndexTabix
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(BedIndexTabix.class);
	
	
	
	private static class BedDataCodec extends  AbstractDataCodec<String>
			{
			@Override
			public String decode(DataInputStream dis)
					throws IOException
				{
				try {
					String s= AbstractDataCodec.readString(dis);
					return s;
				} catch (Exception e) {
					return null;
					}
				}
			@Override
			public void encode(DataOutputStream dos, String s) throws IOException {
				AbstractDataCodec.writeString(dos, s);
					}
			@Override
			public BedDataCodec clone() {
					return new BedDataCodec();
					}
			};

	
	public BedIndexTabix()
		{
		}
	 @Override
		public  Command createCommand() {
				return new MyCommand();
			}
			 
		public static class MyCommand extends AbstractBedIndexTabix.AbstractBedIndexTabixCommand
			 	{
			private int bedLineCount=0;
			private BedLineCodec bedCodec=new BedLineCodec();
	
	
	public int getBedLineCount()
		{
		return bedLineCount;
		}
	
	protected void run(LineIterator in) throws IOException
		{
		this.bedLineCount=0;
		File out=getOutputFile();
		if(out==null)
			{
			throw new IOException(getName()+" output file undefined.");
			}
		if(!out.getName().endsWith(".bed.gz"))
			{
			throw new IOException("output file should en with .bed.gz but got "+out);
			}
		File tbi = new File(getOutputFile().getPath()+TabixUtils.STANDARD_INDEX_EXTENSION);
		BlockCompressedOutputStream writer=null;
		SortingCollection<String> sorter=null;
		final Comparator<String> comparator=new Comparator<String>	()
			{
			@Override
			public int compare(String o1, String o2)
				{
				BedLine bed1 = bedCodec.decode(o1);
				BedLine bed2 = bedCodec.decode(o2);
				int i= bed1.getContig().compareTo(bed2.getContig());
				if(i!=0) return i;
				i= bed1.getStart() - bed2.getStart();
				if(i!=0) return i;
				i= bed1.getEnd() - bed2.getEnd();
				if(i!=0) return i;
				return o1.compareTo(o2);
				}
			};
		CloseableIterator<String> iter=null;

		try {
			TabixIndexCreator indexCreator=new TabixIndexCreator(TabixFormat.BED);
			LOG.info("Opening"+getOutputFile());
			writer = new BlockCompressedOutputStream(getOutputFile());
			
			StringBuilder header=new StringBuilder();
			while(in.hasNext())
				{
				String h=in.peek();
				if(!BedLine.isBedHeader(h)) break;
				header.append(in.next()).append('\n');
				}
			//write header
			if(header.length()>0)
				{
				LOG.info("Writing header");
				writer.write(header.toString().getBytes());
				}
			
			if(this.sort)
				{
				LOG.info("Sorting");
				
				sorter =  SortingCollection.newInstance(
		                        String.class,
		                        new BedDataCodec(),
		                        comparator,
		                        getMaxRecordsInRam(),
		                        getTmpDirectories()
		                        );
				while(in.hasNext())
					{
					String line = in.next();
					BedLine bed = bedCodec.decode(line);
					if(bed==null) continue;
					sorter.add(line);
					}
				sorter.doneAdding();
				sorter.setDestructiveIteration(true);
				iter= sorter.iterator();
				long filePosition= writer.getFilePointer();
				while(iter.hasNext())
					{
					String line = iter.next();
					BedLine bed = this.bedCodec.decode(line);
					writer.write(line.getBytes());
					writer.write('\n');
					this.bedLineCount++;
					indexCreator.addFeature(bed, filePosition);
					filePosition = writer.getFilePointer();
					}
				sorter.cleanup();
				}
			else
				{
				long filePosition= writer.getFilePointer();
				while(in.hasNext())
					{
					String line = in.next();
					BedLine bed = this.bedCodec.decode(line);
					if(bed==null) continue;
					writer.write(line.getBytes());
					writer.write('\n');
					this.bedLineCount++;
					indexCreator.addFeature(bed, filePosition);
					filePosition = writer.getFilePointer();
					}
				}
			writer.flush();
			LOG.info("Creating index");
			Index index = indexCreator.finalizeIndex(writer.getFilePointer());
			LOG.info("Writing index to "+tbi+ " using "+index.getClass());
			index.writeBasedOnFeatureFile(getOutputFile());
			writer.close();
			writer=null;
			LOG.info("Done  N="+bedLineCount);
			} 
		catch (Exception e)
			{
			if(getOutputFile().exists() && getOutputFile().isFile())
				{
				LOG.warn("Deleting "+getOutputFile());
				getOutputFile().delete();
				if(tbi.exists() && tbi.isFile()) tbi.delete();
				}
			throw new IOException(e);
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(sorter);
			CloserUtil.close(writer);
			}
		}
	@Override
	protected Collection<Throwable> call(String inputName)
		{
		if(getOutputFile()==null)
			{
			return wrapException("missing output file");
			}
		if(!getOutputFile().getName().endsWith(".bed.gz"))
			{
			return wrapException("filename should end with bed.gz "+getOutputFile());
			}
		
		LineIterator r = null;
		try {
			if(inputName==null)
				{
				r = IOUtils.openStdinForLineIterator();
				}
			else 
				{
				r= IOUtils.openURIForLineIterator(inputName);
				}
			run(r);
			return RETURN_OK;
		} catch (Exception e)
			{
			return wrapException(e);
			}
		finally
			{
			CloserUtil.close(r);
			}
		}
			 	}
	
		
	public static void main(String[] args)
		{
		new BedIndexTabix().instanceMainWithExit(args);
		}
	}
