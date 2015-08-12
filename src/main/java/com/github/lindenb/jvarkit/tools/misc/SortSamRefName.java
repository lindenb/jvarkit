/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Comparator;

import com.github.lindenb.jvarkit.util.picard.AbstractBamWriterProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class SortSamRefName extends AbstractBamWriterProgram
	{
	private int maxRecordsInRAM=50000;
	
	public void setMaxRecordsInRAM(int maxRecordsInRAM) {
		this.maxRecordsInRAM = maxRecordsInRAM;
	}
	public int getMaxRecordsInRAM() {
		return maxRecordsInRAM;
	}
	
	
	private static class RefNameComparator implements
		Comparator<SAMRecord>
		{
		private SAMRecordQueryNameComparator nameCmp=new SAMRecordQueryNameComparator();
		
		@Override
		public int compare(final SAMRecord o1, final SAMRecord o2)
			{
	        final int refIndex1 = o1.getReferenceIndex();
	        final int refIndex2 = o2.getReferenceIndex();
	        final int cmp = refIndex1 - refIndex2;
	        
	        if (cmp != 0)
	        	{
	        	if (refIndex1 == -1 ) return 1;
	        	if (refIndex2 == -1 ) return -1;
	            return cmp;
	        	}
			return nameCmp.compare(o1, o2);
			}
		
		}
	
	private SortSamRefName()
		{
		
		}
	
	@Override
	public String getProgramDescription()
		{
		return "Sort a BAM of contig and then on name";
		}
	
	@Override
	protected String getOnlineDocUrl()
		{
		return DEFAULT_WIKI_PREFIX+"SortSamRefName";
		}

	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -t (dir) add tmp directory. Optional"); 
		out.println(" -n (int) max records in RAM. Default " +getMaxRecordsInRAM()); 
		super.printOptions(out);
		}

	private int doWork( SamReader in) throws IOException
		{
		SAMFileHeader header= in.getFileHeader();
		
		SAMFileWriter out=null;
		SAMRecordIterator iter=null;
		CloseableIterator<SAMRecord> iter2=null;
		SortingCollection<SAMRecord> sorter=null;
		try
			{
			final BAMRecordCodec bamRecordCodec=new BAMRecordCodec(header);
			final RefNameComparator refNameComparator=new RefNameComparator();
			sorter =SortingCollection.newInstance(
					SAMRecord.class,
					bamRecordCodec,
					refNameComparator,
					getMaxRecordsInRAM(),
					getTmpDirectories()
					);
			sorter.setDestructiveIteration(true);
			
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
			iter = in.iterator();
			while(iter.hasNext())
				{
				sorter.add( progress.watch(iter.next()));
				}
			sorter.doneAdding();
			
			SAMFileHeader header2=header.clone();
			header2.addComment(getProgramName()+" "+getVersion()+" "+getProgramCommandLine());
			header2.setSortOrder(SortOrder.unsorted);
			out = openSAMFileWriter(header2, true);
			iter2 = sorter.iterator();
			while(iter2.hasNext())
				{
				out.addAlignment(iter2.next());
				}
			out.close();out=null;
			sorter.cleanup();
			progress.finish();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			
			CloserUtil.close(iter2);
			CloserUtil.close(iter);
			CloserUtil.close(out);
			}
		}
	
	
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "t:n:"))!=-1)
			{
			switch(c)
				{
				case 't': this.addTmpDirectory(new File(opt.getOptArg()));
				case 'n': this.setMaxRecordsInRAM(Integer.parseInt(opt.getOptArg()));
				default: 
					{
					switch(handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default: break;
						}
					}
				}
			}
		
		SamReader in=null;
		try
			{
			SamReaderFactory srf= SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
			if(opt.getOptInd()==args.length)
				{
				in = srf.open(SamInputResource.of(System.in));
				}
			else if(opt.getOptInd()+1==args.length)
				{
				in = srf.open(SamInputResource.of(args[opt.getOptInd()]));
				}
			else
				{
				error("Illegal number of arguments");
				return -1;
				}
			return doWork(in); 
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			}
		}

	
	public static void main(String[] args) throws IOException
		{
		new SortSamRefName().instanceMainWithExit(args);
		}
		

	}
