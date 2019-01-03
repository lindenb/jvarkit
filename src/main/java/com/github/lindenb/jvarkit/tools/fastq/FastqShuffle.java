/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.fastq;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;

/**

BEGIN_DOC

## Synopsis

```
$ cat f.fq | java -jar dist/fastqshuffle.jar [options] 
$ java -jar dist/fastqshuffle.jar [options] f.fq.gz
$ java -jar dist/fastqshuffle.jar [options] f1.fq.gz f2.fq.gz

```


## Example

```bash
$ $ curl -s "https://raw.githubusercontent.com/bigdatagenomics/adam/fff8ae259e8f6958eefd8de9a3ec39d33392fb21/adam-core/src/test/resources/interleaved_fastq_sample1.fq" |\
java -jar dist/fastqshuffle.jar -i

@H06HDADXX130110:1:2103:11970:57672/1
GGATAGGGTTAGGGTTAGGGTTAGGGCTAGGGATAGGGGTAGGGTTGGGGTTGGTCATCGGGTGTTTCTTTGTGTTTGAGGTTGATTATTGTGATGGTTAAGGTATCTAGGTATTGTAAAAGTTGGCTTTTAACTTAGAAAATTATGTCATTCTGTTCACAAGTGTTTAGATTGGTAGATAGGTACTATGCGATCACTTCCATTGGCTGAGAGTTCGATTGATTATGAGCCACGCTAGTGGTTGAGATCT
+
69+26933-:7;;135,53<>7<692(?2=9:**;<=#####################################################################################################################################################################################################################
@H06HDADXX130110:1:2103:11970:57672/2
AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTATCGTCAAACCTTACCTCCTCCCTAGCCTCCACCCTGACCATGACACCAACCATCAGCCTTATAGAAAACCCCAGAGATGCTCTTATCCTATACCACAATTACCCCATAACGAAAGAAAGGACTGAAAACAAATAAGTAAAATTCGTACAAATTATATCTATGAGTATGTCCCTGAGTGTAGGTGTAGGTGCATCC
+
=>:=>@=?<>>??>;:<?<=;<<?>=;:8;=(5)0-6;1:>?<>##############################################################################################################################################################################################################
@H06JUADXX130110:1:1108:6424:55322/1
AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACTCTAACCCTAACCCTAACCCTAACGGTAACCCTTACCCTTACTGTAACGCTTATCCTAAATCAAATTCTTCCTCTTAAGATCGCTGTTAAAATTAATCCTATTAGAACAGGTCTTCTGGCACCAAGTTATGTCAATATCCCTTACTCTAAACATGCCTTGATCTCTCATGCATCACTTCAGCACAGCTCTTATGGATCTAGGATCCTCAGT
+
=>;=?=@@=?@?@@9>7@=?=;=?@>29?=?;=>@;4@*0878;40'=@;(3399@9>7@:A############################################################################################################################################################################################
@H06JUADXX130110:1:1108:6424:55322/2
AGGGATAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGATAGGGCTAGGGTTAGGGATAGGGATAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTATCGATAGGGATAGGGATAGGGATAGAGTTAGGGCTATGGGTAGGGTTAGAGTCAGGGAAAGAGATAGGGATGGAGATGGGGTTAAAAAGAAGTCAAGGAATTAAGGTAGGGAAACGGTTCGAGATCTGTAAAGGGCAACGA
+
>>;>*9?:@??@@????@????>@?>>@>@?>?????@@???????=<??8;*;:>?;+A?@?>89?@######################################################################################################################################################################################
@H06HDADXX130110:2:2116:3345:91806/1
GTTAGGGTTAGGGTTGGGTTAGGGTTAGGGTTAGGGTTAGGGGTAGGGTTAGGGTTAGGGGTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGGTAGGGCTAGGGTTAAGGGTAGGGTTAGCGAAAGGGCTGGGGTTAGGGGTGCGGGTACGCGTAGCATTAGGGCTAGAAGTAGGATCTGCAGTGCCTGACCGCGTCTGCGCGGCGACTGCCCAAAGCCTGGGGCCGACTCCAGGCTGAAGCTCAT
+
>=<=???>?>???=??>>8<?><=2=<===1194<?;:?>>?#3==>###########################################################################################################################################################################################################
@H06HDADXX130110:2:2116:3345:91806/2
TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTACCCCTAACCCTAACCCTAACCCTAACCCGTACCCTAAACCCAACCCTAACCACAAAGCAAATCCCAACCTTAACCGGAACCCGAAATCTCGCAGCAAATCTGCAGTAGAGACGCAGACTCAACCATGCGTCTATTAGTACGCATTATCATTGCCTCATGCTTCTTAAGTACAGAGAGATGAC
+
==;<?>@@@<>>@??<>>???<=>>?>:><@?4=:>7=5=>:<=@;'@A?########################################################################################################################################################################################################
```


END_DOC
 *
 */
@Program(name="fastqshuffle",
	description="Shuffle Fastq files",
	keywords="fastq"
	)
public class FastqShuffle extends Launcher
	{
	private static final Logger LOG = Logger.build(FastqShuffle.class).make();

	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File fileout = null;

	@Parameter(names={"-i"},description="single input is paired reads interleaved.(optional)")
	private boolean interleaved_input=false;

	@Parameter(names={"-r"},description="random seed ",converter=Launcher.RandomConverter.class)
	private Random random=new Random(-1L);
	
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
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
			final int i= o1.random < o2.random ? -1:  o1.random > o2.random ? 1: 0;
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
		public OneRead decode(final DataInputStream dis) throws IOException
			{
			final OneRead r=new OneRead();
			try {
				r.random = dis.readLong();
			} catch (final IOException e) {
				return null;
				}
			r.index = dis.readLong();
			r.first=readFastqRecord(dis);
			return r;
			}
		@Override
		public void encode(final DataOutputStream dos,final OneRead r)
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
		public TwoReads decode(final DataInputStream dis) throws IOException
			{
			final TwoReads r=new TwoReads();
			try {
				r.random = dis.readLong();
			} catch (final IOException e) {
				return null;
				}
			r.index = dis.readLong();
			r.first=readFastqRecord(dis);
			r.second=readFastqRecord(dis);
			return r;
			}
		@Override
		public void encode(final DataOutputStream dos,final TwoReads r)
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

	
	private static FastqRecord readFastqRecord(final DataInputStream dis)  throws IOException
		{
		final String seqHeader=dis.readUTF();
		final String seqLine=dis.readUTF();
		final String qualHeader=dis.readUTF();
		final String qualLine=dis.readUTF();
		return new FastqRecord(seqHeader, seqLine, qualHeader, qualLine);
		}
	
	private static String notNull(final String s)
		{
		return s==null?"":s;
		}	
	
	private static void writeFastqRecord(final DataOutputStream dos,final FastqRecord r)  throws IOException
		{
		dos.writeUTF(notNull(r.getReadName()));
		dos.writeUTF(notNull(r.getReadString()));
		dos.writeUTF(notNull(r.getBaseQualityHeader()));
		dos.writeUTF(notNull(r.getBaseQualityString()));
		}

	
	public FastqShuffle()
		{
		}
	
	
	private void runPaired(final FastqReader r1, final FastqReader r2,final FastqWriter w1) throws IOException
		{
		long nReads=0;
		final SortingCollection<TwoReads> sorting= SortingCollection.newInstance(
				TwoReads.class,
				new TwoReadsCodec(),
				new TwoReadsCompare(),
				this.writingSortingCollection.getMaxRecordsInRam(),
				this.writingSortingCollection.getTmpPaths()
				);
		sorting.setDestructiveIteration(true);
		while(r1.hasNext())
			{
			final TwoReads p=new TwoReads();
			p.random=this.random.nextLong();
			p.index=nReads;
			p.first=r1.next();
			if(r2!=null)
				{
				if(!r2.hasNext())  throw new IOException("fastq.paired.read.missing");
				p.second=r2.next();
				}
			else
				{
				if(!r1.hasNext())  throw new IOException("fastq.paired.read.missing");
				p.second=r1.next();
				}
			
			if((++nReads)%this.writingSortingCollection.getMaxRecordsInRam()==0)
				{
				LOG.info("Read "+nReads+" reads");
				}

			
			sorting.add(p);
			}
		if(r2!=null && r2.hasNext()) throw new IOException("fastq.paired.read.missing");
		sorting.doneAdding();
		final  CloseableIterator<TwoReads> iter=sorting.iterator();
		
		while(iter.hasNext())
			{
			final TwoReads p=iter.next();
			w1.write(p.first);
			w1.write(p.second);
			}
		
	
		CloserUtil.close(iter);
		sorting.cleanup();
		}
	

	
	private void runSingle(final FastqReader r1,final FastqWriter w1) throws IOException
		{
		long nReads=0;
		final  SortingCollection<OneRead> sorting= SortingCollection.newInstance(
				OneRead.class,
				new OneReadCodec(),
				new OneReadCompare(),
				this.writingSortingCollection.getMaxRecordsInRam(),
				this.writingSortingCollection.getTmpPaths()
				);
		sorting.setDestructiveIteration(true);
		while(r1.hasNext())
			{
			final OneRead r=new OneRead();
			r.random=this.random.nextLong();
			r.index=nReads;
			r.first=r1.next();
			
			if((++nReads)%this.writingSortingCollection.getMaxRecordsInRam()==0)
				{
				LOG.info("Read "+nReads+" reads");
				}

			
			sorting.add(r);
			}
		sorting.doneAdding();
		
		final CloseableIterator<OneRead> iter=sorting.iterator();
		while(iter.hasNext())
			{
			final OneRead p=iter.next();
			w1.write(p.first);
			}
		
	
		CloserUtil.close(iter);
		sorting.cleanup();
		}
	
	@Override
	public int doWork(final List<String> args) {
		FastqReader r1=null;
		FastqReader r2=null;
		FastqWriter w=null;
		
		try
			{
			if(fileout==null)
				{
				w=	new BasicFastqWriter(stdout());
				}
			else
				{
				w=	new BasicFastqWriter(this.fileout);
				}
			
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				r1=new FourLinesFastqReader(stdin());
				if(interleaved_input)
					{
					runPaired(r1, null,w);
					}
				else
					{
					runSingle(r1,w);
					}
				}
			else if(args.size()==1)
				{
				r1=new FourLinesFastqReader(new File(args.get(0)));

				if(interleaved_input)
					{
					runPaired(r1, null,w);
					}
				else
					{
					runSingle(r1,w);
					}
				
				}
			else if(args.size()==2)
				{
				r1=new FourLinesFastqReader(new File(args.get(0)));
				r2=new FourLinesFastqReader(new File(args.get(1)));
				runPaired(r1, r2,w);
				}
			else
				{
				LOG.error("illegal.number.of.arguments");
				return -1;
				}
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
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
