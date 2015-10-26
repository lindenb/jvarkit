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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Comparator;
import java.util.Random;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;



public class VCFShuffle extends AbstractVCFShuffle
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(VcfSetSequenceDictionary.class);

	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractVCFShuffle.AbstractVCFShuffleCommand
		{		
	
	private static class RLine
		{
		long rand;
		String line;
		}
	
	private static class RLineCodec
		extends AbstractDataCodec<RLine>
		{
		@Override
		public RLine decode(DataInputStream dis) throws IOException
			{
			RLine r=new RLine();
			try
				{
				r.rand=dis.readLong();
				}
			catch(IOException err)
				{
				return null;
				}
			r.line=readString(dis);
			return r;
			}
		@Override
		public void encode(DataOutputStream dos, RLine object)
				throws IOException {
			dos.writeLong(object.rand);
			writeString(dos,object.line);
			}
		@Override
		public AbstractDataCodec<RLine> clone() {
			return new RLineCodec();
			}
		}
	
	private static class RLineCmp
		implements Comparator<RLine>
		{
		@Override
		public int compare(RLine o1, RLine o2) {
			int i= o1.rand<o2.rand?-1:o1.rand>o2.rand?1:0;
			if(i!=0) return i;
			return o1.line.compareTo(o2.line);
			}
		}
	
	
	@Override
		protected Collection<Throwable> doVcfToVcf(String inputName)
				throws Exception {
		SortingCollection<RLine> shuffled=null;
		VariantContextWriter out=null;
		BufferedReader lr=null;
		try
			{
			if(inputName==null)
				{
				LOG.info("reading from stdin.");
				lr=new BufferedReader(new InputStreamReader(stdin()));
				}
			else
				{
				lr=IOUtils.openURIForBufferedReading(inputName);
				}
			
			out = super.openVariantContextWriter();
			
			
			Random random=new Random(this.seed);
			
			VCFUtils.CodecAndHeader cah=VCFUtils.parseHeader(lr);
			VCFHeader header=cah.header;
			
			addMetaData(header);
			out.writeHeader(header);
			LOG.info("shuffling");
			
			shuffled=SortingCollection.newInstance(
					RLine.class,
					new RLineCodec(),
					new RLineCmp(),
					getMaxRecordsInRam(),
					getTmpDirectories()
					);
			shuffled.setDestructiveIteration(true);
			String line;
			while((line= lr.readLine())!=null)
				{
				RLine rLine=new RLine();
				rLine.rand=random.nextLong();
				rLine.line=line;
				shuffled.add(rLine);
				}
			shuffled.doneAdding();
			LOG.info("done shuffling");
			
			CloseableIterator<RLine> iter=shuffled.iterator();
			while(iter.hasNext())
				{
				VariantContext ctx=cah.codec.decode(iter.next().line);
				out.add(ctx);
				if(out.checkError()) break;
				}
			LOG.info("Done");
			
			return RETURN_OK;
			}
	catch(Exception err)
		{
		return wrapException(err);
		}
	finally
		{
		if(shuffled!=null) shuffled.cleanup();
		CloserUtil.close(shuffled);
		CloserUtil.close(lr);
		CloserUtil.close(out);
		}
	}
	
	@Override
	public Collection<Throwable> initializeKnime() {
		if(this.seed==-1L) this.seed=System.currentTimeMillis();
		return super.initializeKnime();
		}
	@Override
		protected Collection<Throwable> call(String inputName) throws Exception {
			return doVcfToVcf(inputName);
		}
	

		}
	
	public static void main(String[] args)
		{
		new VCFShuffle().instanceMainWithExit(args);
		}
	}
