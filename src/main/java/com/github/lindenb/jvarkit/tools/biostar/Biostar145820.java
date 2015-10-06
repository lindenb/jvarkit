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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Random;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class Biostar145820 extends AbstractBiostar145820
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Biostar145820.class);

	
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
	public Command createCommand() {
		return new MyCommand();
		}
	static private class MyCommand extends AbstractBiostar145820.AbstractBiostar145820Command
		{
		@Override
		protected Collection<Throwable> call(String inputName) throws Exception
			{
			if(super.count<-1L) // -1 == infinite
				{
				return wrapException("Bad count:"+count);
				}
			SamReader samReader=null;
			SAMRecordIterator iter=null;
			SAMFileWriter samWriter=null;
			Random random=new Random();
			CloseableIterator<RandSamRecord> iter2=null;
			try
				{
				samReader = openSamReader(inputName);
				SAMFileHeader header=samReader.getFileHeader();			
				header=header.clone();
				header.setSortOrder(SortOrder.unsorted);
				header.addComment("Processed with "+getName()+" : "+getProgramCommandLine());
				
				samWriter = super.openSAMFileWriter(header, true);
				
				iter=samReader.iterator();
				SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
				
				SortingCollection<RandSamRecord> sorter=SortingCollection.newInstance(
						RandSamRecord.class,
						new RandSamRecordCodec(header),
						new RandSamRecordComparator(), 
						super.maxRecordsInRam,
						super.tmpDirs
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
				if(super.count==-1)
					{
					while(iter2.hasNext())
						{
						samWriter.addAlignment(iter2.next().samRecord);
						}
					}
				else
					{
					while(iter2.hasNext() && count>0)
						{
						samWriter.addAlignment(iter2.next().samRecord);
						super.count--;
						}
					}
				iter2.close();iter2=null;
				sorter.cleanup();
				progress.finish();
				return Collections.emptyList();
				}
			catch (Exception e)
				{
				LOG.error(e);
				return wrapException(e);
				}
			finally
				{
				CloserUtil.close(iter);
				CloserUtil.close(iter2);
				CloserUtil.close(samReader);
				CloserUtil.close(samWriter);
				}
			}
		}

	
	public static void main(String[] args) throws IOException
		{
		new Biostar145820().instanceMainWithExit(args);
		}
		

	}
