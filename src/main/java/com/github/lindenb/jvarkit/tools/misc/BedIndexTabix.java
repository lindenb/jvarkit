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
import java.io.PrintStream;
import java.util.Comparator;
import java.util.List;

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
import com.github.lindenb.jvarkit.knime.AbstractKnimeApplication;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;


public class BedIndexTabix
	extends AbstractKnimeApplication
	{
	private boolean sort=false;
	private int MAX_RECORDS_IN_RAM=10000;
	private BedLineCodec bedCodec=new BedLineCodec();
	private int bedLineCount=0;
	
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
	public String getProgramDescription() {
		return "Index and sort a Bed on the fly with Tabix.";
		}
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"BedIndexTabix";
		}
	
	public void setSort(boolean sort) {
		this.sort = sort;
		}
	
	public void setMaxRecordsInRam(int n) {
		this.MAX_RECORDS_IN_RAM = n;
		}
	
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
			throw new IOException(getProgramName()+" output file undefined.");
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
			info("Opening"+getOutputFile());
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
				info("Writing header");
				writer.write(header.toString().getBytes());
				}
			
			if(this.sort)
				{
				info("Sorting");
				
				sorter =  SortingCollection.newInstance(
		                        String.class,
		                        new BedDataCodec(),
		                        comparator,
		                        MAX_RECORDS_IN_RAM,
		                        BedIndexTabix.this.getTmpDirectories().get(0)
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
			info("Creating index");
			Index index = indexCreator.finalizeIndex(writer.getFilePointer());
			info("Writing index to "+tbi+ " using "+index.getClass());
			index.writeBasedOnFeatureFile(getOutputFile());
			writer.close();
			writer=null;
			info("Done  N="+bedLineCount);
			} 
		catch (Exception e)
			{
			if(getOutputFile().exists() && getOutputFile().isFile())
				{
				warning("Deleting "+getOutputFile());
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
	public int initializeKnime() {
		return super.initializeKnime();
		}
	@Override
	public void disposeKnime() {
		super.disposeKnime();
		}
	
	@Override
	public int executeKnime(List<String> args) {
		LineIterator r = null;
		try {
			if(args.isEmpty())
				{
				r = IOUtils.openStdinForLineIterator();
				}
			else if(args.size()==1)
				{
				r= IOUtils.openURIForLineIterator(args.get(0));
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			run(r);
			return 0;
		} catch (Exception e)
			{
			error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			}
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -o (file) filename.bed.gz out. Required.");
		out.println(" -s sort BED prior to saving.");
		out.println(" -n (max-record-in-ram). If defined, bed will be sorted using this cache size. Default:"+MAX_RECORDS_IN_RAM);
		out.println(" -T (tmp).add tmp dir. optional)");
		super.printOptions(out);
		}
		
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "o:n:sT:"))!=-1)
			{
			switch(c)
				{
				case 'T': this.addTmpDirectory(new File(opt.getOptArg()));break;
				case 's': this.setSort(true);break;
				case 'n': this.setMaxRecordsInRam(Integer.parseInt(opt.getOptArg()));break;
				case 'o': this.setOutputFile(new File(opt.getOptArg()));break;
				default: 
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		return mainWork(opt.getOptInd(), args);
		}
		
		
		
	public static void main(String[] args)
		{
		new BedIndexTabix().instanceMainWithExit(args);
		}
	}
