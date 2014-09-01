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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.fastq;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.Random;

import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;

/**
 * FastqShuffle
 * @author lindenb
 *
 */
public class FastqShuffle extends AbstractCommandLineProgram
	{
	private Random random=new Random();
	private int maxRecordsInRAM=500000;
	
	private static class OneRead
		{
		long random;
		long index;
		FastqRecord first;
		}

	
	private static class TwoReads extends OneRead
		{
		FastqRecord second;
		}
	
	private static class OneReadCompare
		implements Comparator<OneRead>
		{
		@Override
		public int compare(final OneRead o1, final OneRead o2) {
			int i= o1.random < o2.random ? -1:  o1.random > o2.random ? 1: 0;
			if(i!=0) return i;
			return  o1.index < o2.index ? -1:  o1.index > o2.index ? 1: 0;
			}
		}
	private static class TwoReadsCompare
	implements Comparator<TwoReads>
		{
		@Override
		public int compare(final TwoReads o1, final TwoReads o2) {
			int i= o1.random < o2.random ? -1:  o1.random > o2.random ? 1: 0;
			if(i!=0) return i;
			return  o1.index < o2.index ? -1:  o1.index > o2.index ? 1: 0;
			}
		}

	
	private static class OneReadCodec extends AbstractDataCodec<OneRead>
		{
		@Override
		public OneRead decode(DataInputStream dis) throws IOException
			{
			OneRead r=new OneRead();
			try {
				r.random = dis.readLong();
			} catch (IOException e) {
				return null;
				}
			r.index = dis.readLong();
			r.first=readFastqRecord(dis);
			return r;
			}
		@Override
		public void encode(DataOutputStream dos, OneRead r)
				throws IOException {
			dos.writeLong(r.random);
			dos.writeLong(r.index);
			writeFastqRecord(dos,r.first);
			}
		@Override
		public AbstractDataCodec<OneRead> clone() {
			return new OneReadCodec();
			}
		}
	
	private static class TwoReadsCodec extends AbstractDataCodec<TwoReads>
		{
		@Override
		public TwoReads decode(DataInputStream dis) throws IOException
			{
			TwoReads r=new TwoReads();
			try {
				r.random = dis.readLong();
			} catch (IOException e) {
				return null;
				}
			r.index = dis.readLong();
			r.first=readFastqRecord(dis);
			r.second=readFastqRecord(dis);
			return r;
			}
		@Override
		public void encode(DataOutputStream dos, TwoReads r)
				throws IOException {
			dos.writeLong(r.random);
			dos.writeLong(r.index);
			writeFastqRecord(dos,r.first);
			writeFastqRecord(dos,r.second);
			}
		@Override
		public AbstractDataCodec<TwoReads> clone() {
			return new TwoReadsCodec();
			}
		}

	
	private static FastqRecord readFastqRecord(DataInputStream dis)  throws IOException
		{
		String seqHeader=dis.readUTF();
		String seqLine=dis.readUTF();
		String qualHeader=dis.readUTF();
		String qualLine=dis.readUTF();
		return new FastqRecord(seqHeader, seqLine, qualHeader, qualLine);
		}
	
	private static String notNull(String s)
		{
		return s==null?"":s;
		}	
	
	private static void writeFastqRecord(DataOutputStream dos,final FastqRecord r)  throws IOException
		{
		dos.writeUTF(notNull(r.getReadHeader()));
		dos.writeUTF(notNull(r.getReadString()));
		dos.writeUTF(notNull(r.getBaseQualityHeader()));
		dos.writeUTF(notNull(r.getBaseQualityString()));
		}

	
	private FastqShuffle()
		{
		
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/FastqShuffle";
		}
	
	@Override
	public String getProgramDescription() {
		return "Shuffle Fastq files.";
		}
	
	
	private void runPaired(FastqReader r1, FastqReader r2,FastqWriter w1) throws IOException
		{
		long nReads=0;
		SortingCollection<TwoReads> sorting= SortingCollection.newInstance(
				TwoReads.class,
				new TwoReadsCodec(),
				new TwoReadsCompare(),
				maxRecordsInRAM,
				getTmpDirectories()
				);
		sorting.setDestructiveIteration(true);
		while(r1.hasNext())
			{
			if(r2!=null)
				{
				if(!r2.hasNext()) throw new IOException(getMessageBundle("fastq.paired.read.missing"));
				}
			else
				{
				if(!r1.hasNext()) throw new IOException(getMessageBundle("fastq.paired.read.missing"));
				}
			TwoReads p=new TwoReads();
			p.random=this.random.nextLong();
			p.index=nReads;
			p.first=r1.next();
			p.second=(r2==null?r1.next():r2.next());
			
			
			if((++nReads)%this.maxRecordsInRAM==0)
				{
				info("Read "+nReads+" reads");
				}

			
			sorting.add(p);
			}
		if(r2!=null && r2.hasNext()) throw new IOException(getMessageBundle("fastq.paired.read.missing"));
		sorting.doneAdding();
		CloseableIterator<TwoReads> iter=sorting.iterator();
		
		while(iter.hasNext())
			{
			TwoReads p=iter.next();
			w1.write(p.first);
			w1.write(p.second);
			}
		
	
		CloserUtil.close(iter);
		sorting.cleanup();
		}
	

	
	private void runSingle(FastqReader r1,FastqWriter w1) throws IOException
		{
		long nReads=0;
		SortingCollection<OneRead> sorting= SortingCollection.newInstance(
				OneRead.class,
				new OneReadCodec(),
				new OneReadCompare(),
				maxRecordsInRAM,
				getTmpDirectories()
				);
		sorting.setDestructiveIteration(true);
		while(r1.hasNext())
			{
			OneRead r=new OneRead();
			r.random=this.random.nextLong();
			r.index=nReads;
			r.first=r1.next();
			
			if((++nReads)%this.maxRecordsInRAM==0)
				{
				info("Read "+nReads+" reads");
				}

			
			sorting.add(r);
			}
		CloseableIterator<OneRead> iter=sorting.iterator();
		while(iter.hasNext())
			{
			OneRead p=iter.next();
			w1.write(p.first);
			}
		
	
		CloserUtil.close(iter);
		sorting.cleanup();
		}
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -i single input is paired reads interleaved.(optional).");
		out.println(" -o fastq interleaved output. (optional).");
		out.println(" -N (int) "+getMessageBundle("max.records.in.ram")+" default:"+maxRecordsInRAM+". Optional.");
		out.println(" -T (dir) "+getMessageBundle("add.tmp.dir")+" Optional.");
		super.printOptions(out);
		}

	@Override
	public int doWork(String[] args)
		{
		File fileout=null;
		boolean interleaved_input=false;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:N:T:i"))!=-1)
			{
			switch(c)
				{
				case 'i': interleaved_input=true;break;
				case 'o': fileout=new File(opt.getOptArg());break;
				case 'N':this.maxRecordsInRAM=Integer.parseInt(opt.getOptArg());break;
				case 'T':this.addTmpDirectory(new File(opt.getOptArg()));break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		FastqReader r1=null;
		FastqReader r2=null;
		FastqWriter w=null;
		
		try
			{
			if(fileout==null)
				{
				w=	new BasicFastqWriter(System.out);
				}
			else
				{
				w=	new BasicFastqWriter(fileout);
				if(fileout.getParentFile()!=null && fileout.getParentFile().exists())
					{
					this.addTmpDirectory( fileout.getParentFile());
					}
				}
			
			if(opt.getOptInd() == args.length)
				{
				info("Reading from stdin");
				r1=new FourLinesFastqReader(System.in);
				if(interleaved_input)
					{
					runPaired(r1, null,w);
					}
				else
					{
					runSingle(r1,w);
					}
				}
			else if(opt.getOptInd()+1==args.length)
				{
				r1=new FourLinesFastqReader(new File(args[opt.getOptInd()]));

				if(interleaved_input)
					{
					runPaired(r1, null,w);
					}
				else
					{
					runSingle(r1,w);
					}
				
				}
			else if(opt.getOptInd()+2==args.length)
				{
				r1=new FourLinesFastqReader(new File(args[opt.getOptInd()  ]));
				r2=new FourLinesFastqReader(new File(args[opt.getOptInd()+1]));
				runPaired(r1, r2,w);
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r1);
			CloserUtil.close(r2);
			CloserUtil.close(w);
			}
		}
	
	public static void main(String[] args)
		{
		new FastqShuffle().instanceMainWithExit(args);
		}

}
