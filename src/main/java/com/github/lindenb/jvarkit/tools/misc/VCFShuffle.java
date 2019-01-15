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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Random;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
/**
BEGIN_DOC

## Example

```
$ java -jar dist/vcfshuffle.jar input.vcf
```


END_DOC
 */
@Program(
	name="vcfshuffle",
	description="Shuffle a VCF",
	keywords={"vcf"}
	)
public class VCFShuffle extends Launcher
	{
	private static final Logger LOG = Logger.build(VCFShuffle.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	
	@Parameter(names={"-N","--seed"},description="random seed. Optional. -1 = time.")
	private long seed = -1L ;

	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
	
	private static class RLine
		{
		long rand;
		String line;
		}
	
	private static class RLineCodec
		extends AbstractDataCodec<RLine>
		{
		@Override
		public RLine decode(final DataInputStream dis) throws IOException
			{
			final RLine r=new RLine();
			try
				{
				r.rand=dis.readLong();
				}
			catch(final EOFException err)
				{
				return null;
				}
			r.line=readString(dis);
			return r;
			}
		@Override
		public void encode(final DataOutputStream dos,final  RLine object)
				throws IOException {
			dos.writeLong(object.rand);
			writeString(dos,object.line);
			}
		@Override
		public AbstractDataCodec<RLine> clone() {
			return new RLineCodec();
			}
		}
		
	public VCFShuffle()
		{
		}
	
	@Override
	public int doWork(final List<String> args) {
		if(seed==-1L) seed= System.currentTimeMillis();
		SortingCollection<RLine> shuffled=null;
		VariantContextWriter out=null;
		BufferedReader lr=null;
		try
			{
			lr = super.openBufferedReader(oneFileOrNull(args));
			out = super.openVariantContextWriter(this.outputFile);
			
			
			final Random random=new Random(this.seed);

			final VCFUtils.CodecAndHeader cah=VCFUtils.parseHeader(lr);
			final VCFHeader header=cah.header;
			super.addMetaData(header);
			out.writeHeader(header);
			LOG.info("shuffling");
			
			shuffled=SortingCollection.newInstance(
					RLine.class,
					new RLineCodec(),
					(o1,o2)->{
						final int i= Long.compare(o1.rand, o2.rand);
						if(i!=0) return i;
						return o1.line.compareTo(o2.line);
						},
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			shuffled.setDestructiveIteration(true);
			String line;
			while((line= lr.readLine())!=null)
				{
				final RLine rLine=new RLine();
				rLine.rand=random.nextLong();
				rLine.line=line;
				shuffled.add(rLine);
				}
			shuffled.doneAdding();
			LOG.info("done shuffling");
			
			final CloseableIterator<RLine> iter=shuffled.iterator();
			while(iter.hasNext())
				{
				final VariantContext ctx=cah.codec.decode(iter.next().line);
				out.add(ctx);
				if(out.checkError()) break;
				}
			return RETURN_OK;
			}
	catch(final Exception err)
		{
		LOG.error(err);
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
	

	public static void main(final String[] args)
		{
		new VCFShuffle().instanceMainWithExit(args);
		}
	}
