package com.github.lindenb.jvarkit.tools.misc;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

/**
BEGIN_DOC


END_DOC
 */
@Program(name="coveragenormalizer",description="normalize BAM coverage" )
public class CoverageNormalizer extends Launcher
	{
	private static Logger LOG=Logger.build(CoverageNormalizer.class).make();

	@Parameter(names="-w",description=" (size) sliding window size")
	private int window_size=100;
	@Parameter(names="-s",description=" (size) sliding window shift")
	private int window_shift=50;
	
	private int min_coverage=10;
	
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();

	
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
	
	private CoverageNormalizer()
		{
		
		}
	
	private int run(SamReader sfr) throws IOException
		{
		SAMSequenceDictionary dictionary=sfr.getFileHeader().getSequenceDictionary();
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dictionary);
		File tmpFile1=File.createTempFile("cov_", ".dat.gz", this.writingSortingCollection.getTmpDirectories().get(0));
		tmpFile1.deleteOnExit();

		LOG.info("Opening tmp File "+tmpFile1);
		GZIPOutputStream gos=null;
		DataInputStream dis=null;
		SortingCollection<Float> median=null;
		try
			{
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
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
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
						System.out.print(dictionary.getSequence(chrom_id).getSequenceName());
						System.out.print('\t');
						System.out.print(i*window_shift);
						System.out.print('\t');
						System.out.print((j-1)*window_shift+window_size);
						System.out.print('\t');
						System.out.print((float)((value_start-minCov)/(double)(maxCov-minCov)));
						System.out.print('\t');
						System.out.print(median_value);
						System.out.println();
						if(System.out.checkError()) break;
						}
					
					i=j;
					if(value_end==null) break;
					value_start=value_end;
					}
			
			 	}
			 CloserUtil.close(dis);

			 
			 
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(gos);
			CloserUtil.close(dis);
			if(tmpFile1!=null) tmpFile1.delete();
			if(median!=null) median.cleanup();
			}
		
		
		}
	
	@Override
	public int doWork(List<String> args)
		{
		SamReader sfr=null;
		try
			{
			sfr = super.openSamReader(oneFileOrNull(args));
			
			if(sfr.getFileHeader().getSortOrder()!=SortOrder.coordinate)
				{
				LOG.error("BAM file is not sorted.");
				return -1;
				}
			
			run(sfr);
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(sfr);
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
