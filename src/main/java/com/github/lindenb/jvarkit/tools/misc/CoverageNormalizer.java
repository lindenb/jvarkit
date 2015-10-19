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




*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class CoverageNormalizer extends AbstractCoverageNormalizer
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(ConvertBedChromosomes.class);

	
	private static class FloatCodec extends AbstractDataCodec<Float>
		{
		@Override
		public Float decode(DataInputStream dis) throws IOException
			{
			try { return dis.readFloat();}
			catch(Exception err) { return null;}
			}
		@Override
		public void encode(DataOutputStream dos, Float v)
				throws IOException {
			dos.writeFloat(v);
			}
		@Override
		public AbstractDataCodec<Float> clone() {
			return new FloatCodec();
			}
		}
	private static class FloatCmp
		implements Comparator<Float>
		{
		@Override
		public int compare(Float f1, Float f2) {
			return f1.compareTo(f2);
			}
		}

	 @Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  class MyCommand extends AbstractCoverageNormalizer.AbstractCoverageNormalizerCommand
	 	{		
		private int min_coverage=10;
		
		@Override
		protected Collection<Throwable> call(String inputName)
			{
			
			File tmpFile1=null;
	
			LOG.info("Opening tmp File "+tmpFile1);
			GZIPOutputStream gos=null;
			DataInputStream dis=null;
			SortingCollection<Float> median=null;
			SamReader sfr = null;
			PrintStream out = null;
			try
				{
				out  = openFileOrStdoutAsPrintStream();
				sfr = openSamReader(inputName);
				SAMSequenceDictionary dictionary=sfr.getFileHeader().getSequenceDictionary();
				SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dictionary);
				tmpFile1 = File.createTempFile("cov_", ".dat.gz", getTmpDirectories().get(0));
				tmpFile1.deleteOnExit();
				
				gos=new GZIPOutputStream(new FileOutputStream(tmpFile1));
				DataOutputStream daos=new DataOutputStream(gos);
				SAMRecordIterator iter=sfr.iterator();
				int curr_tid=-1;
				short array[]=null;
				float minCov=Float.MAX_VALUE;
				float maxCov=0;
				int num_written[]=new int[dictionary.size()];
				Arrays.fill(num_written, 0);
				
				for(;;)
					{
					SAMRecord rec=null;
					if(iter.hasNext())
						{
						rec=iter.next();
						if(rec.getReadUnmappedFlag()) continue;
						progress.watch(rec);
						if(rec.isSecondaryOrSupplementary()) continue;
						}
					if(rec==null || curr_tid!=rec.getReferenceIndex())
						{
						if(curr_tid!=-1)
							{
							LOG.info("Writing data for chromosome "+dictionary.getSequence(curr_tid).getSequenceName());
							for(int i=0;
									i + window_size <= array.length;
									i+=window_shift
									)
								{
								float sum=0;
								for(int j=0;j< window_size;++j)
									{
									sum+=array[i+j];
									}
								float v=sum/window_size;
								daos.writeFloat(v);
								minCov=(float)Math.min(minCov, v);
								maxCov=(float)Math.max(maxCov, v);
								num_written[curr_tid]++;
								}
							LOG.info("End writing data N="+num_written[curr_tid]);
							}
						if(rec==null) break;
						curr_tid=rec.getReferenceIndex();
						SAMSequenceRecord ssr=dictionary.getSequence(curr_tid);
						LOG.info("allocate "+ssr.getSequenceLength()+" for "+ssr.getSequenceName());
						array=null;
						System.gc();
						array=new short[ssr.getSequenceLength()];
						LOG.info("done: allocate.");
						Arrays.fill(array, (short)0);
						}
					for(int i=rec.getAlignmentStart();i< rec.getAlignmentEnd() && i< array.length ; ++i)
						{
						array[i]=(short)Math.min((int)Short.MAX_VALUE,1+(int)array[i]);
						}
					}
				array=null;
				LOG.info("Closing BAM");
				CloserUtil.close(sfr);
				
				LOG.info("Closing "+tmpFile1);
				daos.flush();
				gos.finish();
				gos.flush();
				gos.close();
				
				//start normalizing min/max find median value
				long nWritten=0L;
				median=SortingCollection.newInstance(
						Float.class,
						new FloatCodec(),
						new FloatCmp(),
						this.getMaxRecordsInRam(),
						this.getTmpDirectories()
						);
				 median.setDestructiveIteration(true);
				 dis=new DataInputStream(new GZIPInputStream(new FileInputStream(tmpFile1)));
				 for(int n_items:num_written)
				 	{
					for(int i=0;i< n_items;++i)
						{
						float v=dis.readFloat();
						if(v<min_coverage) continue;
						v=(float)((v-minCov)/(double)(maxCov-minCov));
						median.add(v);
						++nWritten;
						}
				 	}
				 median.doneAdding();
				 CloserUtil.close(dis);
	
				 //get median
				 float median_value=0f;
				 long half=nWritten/2L;
				 CloseableIterator<Float> iterFloat=median.iterator();
				 while(iterFloat.hasNext() && half>0)
				 	{
					 median_value=iterFloat.next();
					--half;
					if(half<=0)
						{
						LOG.info("median = "+median_value);
						break;
						}
				 	}
				 CloserUtil.close(iterFloat);
				 median.cleanup();
				 median=null;
				
				 progress=new SAMSequenceDictionaryProgress(dictionary);
				 //dump data
				 dis=new DataInputStream(new GZIPInputStream(new FileInputStream(tmpFile1)));
				 for(int chrom_id=0;
						 chrom_id< num_written.length;
						 ++chrom_id
						 )
				 	{
					int n_items=num_written[chrom_id];
					int i=0;
					Float value_start=null;
					while(i< n_items)
						{
						if(value_start==null)
							{
							value_start=dis.readFloat();
							}
						int j=i+1;
						Float value_end=null;
						while(j < n_items)
							{
							value_end=dis.readFloat();
							if(value_start.intValue()==value_end.intValue())
								{
								++j;
								}
							else
								{
								break;
								}
							}
						progress.watch(chrom_id, i*window_shift);
						if(value_start>=min_coverage)
							{
							out.print(dictionary.getSequence(chrom_id).getSequenceName());
							out.print('\t');
							out.print(i*window_shift);
							out.print('\t');
							out.print((j-1)*window_shift+window_size);
							out.print('\t');
							out.print((float)((value_start-minCov)/(double)(maxCov-minCov)));
							out.print('\t');
							out.print(median_value);
							out.println();
							if(out.checkError()) break;
							}
						
						i=j;
						if(value_end==null) break;
						value_start=value_end;
						}
				
				 	}
				 CloserUtil.close(dis);
				 
				return RETURN_OK;
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(gos);
				CloserUtil.close(dis);
				CloserUtil.close(sfr);
				CloserUtil.close(out);
				if(tmpFile1!=null) tmpFile1.delete();
				if(median!=null) median.cleanup();
				}			
			}
	 	}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new CoverageNormalizer().instanceMainWithExit(args);
		}		

}
