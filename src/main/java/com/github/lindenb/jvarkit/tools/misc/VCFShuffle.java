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
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.knime.AbstractKnimeApplication;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;



public class VCFShuffle extends AbstractKnimeApplication
	{
	private int maxRecordsInRAM=50000;
	private long seed=System.currentTimeMillis();
	private int variantCount=0;
	
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

	
	public VCFShuffle()
		{
		}
	

	@Override
	public String getProgramDescription() {
		return "Shuffle a VCF";
		}
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"VcfShuffle";
		}
	
	@Override
	public int initializeKnime() {
		return super.initializeKnime();
		}

	@Override
	public void disposeKnime() {
		super.disposeKnime();
		}
	
	public int getVariantCount() {
		return variantCount;
	}
	
	@Override
	public int executeKnime(List<String> args)
		{
		SortingCollection<RLine> shuffled=null;
		VariantContextWriter out=null;
		BufferedReader lr=null;
		try
			{
			if(args.isEmpty())
				{
				info("reading from stdin.");
				lr=IOUtils.openStdinForBufferedReader();
				}
			else if(args.size()==1)
				{
				String filename=args.get(0);
				info("reading from "+filename);
				lr=IOUtils.openURIForBufferedReading(filename);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			if(getOutputFile()==null)
				{
				out = VCFUtils.createVariantContextWriterToStdout();
				}
			else
				{
				out = VCFUtils.createVariantContextWriter(getOutputFile());
				}
			
			File tmpFile = super.getTmpDirectories().get(0);
			Random random=new Random(this.seed);
			
			VCFUtils.CodecAndHeader cah=VCFUtils.parseHeader(lr);
			VCFHeader header=cah.header;
			
			header.addMetaDataLine(new VCFHeaderLine("VCFShuffle.Version",String.valueOf(getVersion())));
			out.writeHeader(header);
			info("shuffling");
			
			shuffled=SortingCollection.newInstance(
					RLine.class,
					new RLineCodec(),
					new RLineCmp(),
					maxRecordsInRAM,
					tmpFile
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
			info("done shuffling");
			
			this.variantCount=0;
			CloseableIterator<RLine> iter=shuffled.iterator();
			while(iter.hasNext())
				{
				VariantContext ctx=cah.codec.decode(iter.next().line);
				out.add(ctx);
				++variantCount;
				if(checkOutputError()) break;
				}
			info("Done N="+variantCount);
			
			return 0;
			}
	catch(Exception err)
		{
		error(err);
		return -1;
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
	public void printOptions(PrintStream out) {
		out.println(" -T (dir) tmp directory. Optional.");
		out.println(" -N (long) random seed. Optional.");
		out.println(" -m (int) max records in ram. Optional");
		out.println(" -o (file) output file (default stdout)");
		super.printOptions(out);
		}
	
	public void setSeed(long seed) {
		this.seed = seed;
	}
	
	public void setMaxRecordsInRAM(int maxRecordsInRAM) {
		this.maxRecordsInRAM = maxRecordsInRAM;
	}
	

	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "T:N:m:o:"))!=-1)
			{
			switch(c)
				{
				case 'o': setOutputFile(opt.getOptArg()); break;
				case 'T': addTmpDirectory(new File(opt.getOptArg())); break;
				case 'N': setSeed(Long.parseLong(opt.getOptArg())); break;
				case 'm': setMaxRecordsInRAM(Math.max(10, Integer.parseInt(opt.getOptArg()))); break;
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
		new VCFShuffle().instanceMainWithExit(args);
		}
	}
