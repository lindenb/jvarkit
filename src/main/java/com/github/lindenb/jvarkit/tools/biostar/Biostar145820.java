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
package com.github.lindenb.jvarkit.tools.biostar;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Comparator;
import java.util.Random;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;

public class Biostar145820 extends AbstractCommandLineProgram
	{
	private static class RandSamRecord
		{
		int rand_index;
		SAMRecord samRecord;
		}
	private static class RandSamRecordComparator
		implements Comparator<RandSamRecord>
		{
		SAMRecordQueryNameComparator secondCompare =new SAMRecordQueryNameComparator();
		@Override
		public int compare(RandSamRecord o1, RandSamRecord o2)
			{
			long i = (long)o1.rand_index - (long)o2.rand_index;
			if(i!=0L) return (i<0?-1:1);
			return secondCompare.compare(o1.samRecord, o2.samRecord);
			}
		}
	
	private static class RandSamRecordCodec
		implements SortingCollection.Codec<RandSamRecord>
		{
	    private final BinaryCodec binaryCodec = new BinaryCodec();
	    private BAMRecordCodec bamRecordCodec;
	    private SAMFileHeader header;
		RandSamRecordCodec(SAMFileHeader header)
			{
			this.header=header;
			this.bamRecordCodec=new BAMRecordCodec(this.header);
			}
		
		@Override
		public RandSamRecord decode()
			{
			int r;
			  try {
		             r = this.binaryCodec.readInt();
		        }
		     catch (Exception e) {
		            return null;
		        }
			  RandSamRecord o =new RandSamRecord();
			  o.rand_index=r;
			  o.samRecord = this.bamRecordCodec.decode();
			  return o;
			}
		
		@Override
		public void encode(RandSamRecord val) {
			this.binaryCodec.writeInt(val.rand_index);
			this.bamRecordCodec.encode(val.samRecord);
			}
		@Override
		public void setInputStream(InputStream is) {
			this.binaryCodec.setInputStream(is);
			this.bamRecordCodec.setInputStream(is);
			}
		@Override
		public void setOutputStream(OutputStream os) {
			this.binaryCodec.setOutputStream(os);
			this.bamRecordCodec.setOutputStream(os);
			}
		
		@Override
		public RandSamRecordCodec clone()  {
			return new RandSamRecordCodec(this.header);
			}
		}
	
	
	@Override
	public String getProgramDescription()
		{
		return "subsample BAM to fixed number of alignments. see https://www.biostars.org/p/145820/";
		}
	
	@Override
	protected String getOnlineDocUrl()
		{
		return DEFAULT_WIKI_PREFIX+"Biostar145820";
		}

	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -n (int) number of reads"); 
		out.println(" -o (file) output file (default stdout)"); 
		out.println(" -N (int) max records in ram (optional)"); 
		out.println(" -b force binary for stdout (optional)"); 
		out.println(" -T (dir) add tmp directory (optional)"); 
		super.printOptions(out);
		}

	@Override
	public int doWork(String[] args)
		{
		boolean compressed=false;
		int maxRecordsInRAM=100000;
		long count=10L;
		File fileout=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "o:n:N:T:b"))!=-1)
			{
			switch(c)
				{
				case 'b': compressed=true;break;
				case 'N': maxRecordsInRAM = Integer.parseInt(opt.getOptArg());break;				
				case 'n': count = Long.parseLong(opt.getOptArg());break;				
				case 'o': fileout = new File(opt.getOptArg());break;				
				case 'T': this.addTmpDirectory(new File(opt.getOptArg()));break;				
				default: 
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default: break;
						}
					}
				}
			}
		if(count<=0L)
			{
			error("Bad count");
			return -1;
			}
		SamReader samReader=null;
		SAMRecordIterator iter=null;
		SAMFileWriter samWriter=null;
		Random random=new Random();
		CloseableIterator<RandSamRecord> iter2=null;
		try
			{
			SamFileReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				samReader=SamFileReaderFactory.mewInstance().openStdin();
				}
			else if(opt.getOptInd()+1==args.length)
				{
				File filename=new File(args[opt.getOptInd()]);
				info("Reading from "+filename);
				samReader=SamFileReaderFactory.mewInstance().open(filename);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			SAMFileHeader header=samReader.getFileHeader();
						
			header=header.clone();
			header.setSortOrder(SortOrder.unsorted);
			SAMFileWriterFactory sfw=new SAMFileWriterFactory();
			sfw.setCreateIndex(false);
			sfw.setCreateMd5File(false);
			if(fileout==null)
				{
				if(compressed)
					{
					samWriter=sfw.makeBAMWriter(header,true, System.out);
					}
				else
					{
					samWriter=sfw.makeSAMWriter(header,true, System.out);
					}
				}
			else 
				{
				samWriter=sfw.makeSAMOrBAMWriter(header,true,fileout);
				this.addTmpDirectory(fileout);
				}
			iter=samReader.iterator();
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(samReader.getFileHeader().getSequenceDictionary());
			
			SortingCollection<RandSamRecord> sorter=SortingCollection.newInstance(
					RandSamRecord.class,
					new RandSamRecordCodec(header),
					new RandSamRecordComparator(), 
					maxRecordsInRAM,
					getTmpDirectories()
					);
			sorter.setDestructiveIteration(true);
			while(iter.hasNext())
				{
				RandSamRecord r=new RandSamRecord();
				r.rand_index  = random.nextInt();
				r.samRecord =  progress.watch(iter.next());

				sorter.add(r);
				}
			iter.close();iter=null;
			
			sorter.doneAdding();
			iter2=sorter.iterator();
			while(iter2.hasNext() && count>0)
				{
				samWriter.addAlignment(iter2.next().samRecord);
				count--;
				}
			iter2.close();iter2=null;
			sorter.cleanup();
			progress.finish();
			}
		catch (Exception e)
			{
			error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(iter2);
			CloserUtil.close(samReader);
			CloserUtil.close(samWriter);
			}
		return 0;
		}

	
	public static void main(String[] args) throws IOException
		{
		new Biostar145820().instanceMainWithExit(args);
		}
		

	}
