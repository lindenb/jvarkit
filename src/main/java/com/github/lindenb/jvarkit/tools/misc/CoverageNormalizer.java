package com.github.lindenb.jvarkit.tools.misc;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.SortingCollection;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class CoverageNormalizer extends AbstractCommandLineProgram
	{
	private int maxRecordsInRAM=500000;
	private int window_size=100;
	private int window_shift=50;
	private File tmpDir=null;
	private int min_coverage=10;
	

	
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
	@Override
	public String getProgramDescription() {
		return "normalize BAM coverage";
		}
	
	private int run(SAMFileReader sfr) throws IOException
		{
		sfr.setValidationStringency(ValidationStringency.LENIENT);
		SAMSequenceDictionary dictionary=sfr.getFileHeader().getSequenceDictionary();
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dictionary);
		File tmpFile1=File.createTempFile("cov_", ".dat.gz", this.tmpDir);
		tmpFile1.deleteOnExit();

		info("Opening tmp File "+tmpFile1);
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
						info("Writing data for chromosome "+dictionary.getSequence(curr_tid).getSequenceName());
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
						info("End writing data N="+num_written[curr_tid]);
						}
					if(rec==null) break;
					curr_tid=rec.getReferenceIndex();
					SAMSequenceRecord ssr=dictionary.getSequence(curr_tid);
					info("allocate "+ssr.getSequenceLength()+" for "+ssr.getSequenceName());
					array=null;
					System.gc();
					array=new short[ssr.getSequenceLength()];
					info("done: allocate.");
					Arrays.fill(array, (short)0);
					}
				for(int i=rec.getAlignmentStart();i< rec.getAlignmentEnd() && i< array.length ; ++i)
					{
					array[i]=(short)Math.min((int)Short.MAX_VALUE,1+(int)array[i]);
					}
				}
			array=null;
			info("Closing BAM");
			CloserUtil.close(sfr);
			
			info("Closing "+tmpFile1);
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
					this.maxRecordsInRAM,
					this.tmpDir
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
					info("median = "+median_value);
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
			error(err);
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
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -T (dir) tmpDir");
		out.println(" -w (size) sliding window size");
		out.println(" -s (size) sliding window shift");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"T:w:s:"))!=-1)
			{
			switch(c)
				{
				case 'T':tmpDir=new File(opt.getOptArg());break;
				case 'w':window_size=Math.max(1, Integer.parseInt(opt.getOptArg()));break;
				case 's':window_shift=Math.max(1, Integer.parseInt(opt.getOptArg()));break;
				default:
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		SAMFileReader sfr=null;
		try
			{
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				sfr=new SAMFileReader(System.in);
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				info("Reading from "+filename);
				File f=new File(filename);
				if(tmpDir==null) tmpDir=f.getParentFile();
				sfr=new SAMFileReader(f);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			if(sfr.getFileHeader().getSortOrder()!=SortOrder.coordinate)
				{
				error("BAM file is not sorted.");
				return -1;
				}
			if(tmpDir==null)
				{
				tmpDir=new File(System.getProperty("java.io.tmpdir"));
				info("Setting tmpDir to "+tmpDir);
				}
			run(sfr);
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			if(sfr!=null) sfr.close();
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
