/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.util.Random;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFIterator;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
/**
BEGIN_DOC

## Example

```
$ java -jar dist/vcfshuffle.jar input.vcf
```

## native alternative

```
bcftools view --header-only in.vcf > tmp1.vcf
bcftools view --no-header in.vcf |\
	awk '{printf("%d\t%s\n",int(rand()*10000),$0);}' |\
	sort -t $'\t' -k1,1n -T . |\
	cut -f 1 > tmp2.vcf
	
cat tmp1.vcf tmp2.vcf > shuffled.vcf
```

END_DOC
 */
@Program(
	name="vcfshuffle",
	description="Shuffle a VCF",
	keywords={"vcf"},
	creationDate="20131210",
	modificationDate="20200818"
	)
public class VCFShuffle extends OnePassVcfLauncher
	{
	private static final Logger LOG = Logger.build(VCFShuffle.class).make();
	
	@Parameter(names={"-N","--seed"},description="random seed. Optional. -1 = use current time.")
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
	protected int doVcfToVcf(final String inputName, final VCFIterator in, final VariantContextWriter out) {
		SortingCollection<RLine> shuffled=null;

		try {			
			final Random random=new Random(this.seed);
			final VCFHeader header = in.getHeader();
			final VCFEncoder vcfEncoder = new VCFEncoder(header, false, false);
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

			while(in.hasNext())
				{
				final RLine rLine=new RLine();
				rLine.rand=random.nextLong();
				rLine.line=vcfEncoder.encode(in.next());
				shuffled.add(rLine);
				}
			shuffled.doneAdding();
			
			JVarkitVersion.getInstance().addMetaData(this, header);
			out.writeHeader(header);
			final VCFCodec vcfCodec = new VCFCodec();
			vcfCodec.setVCFHeader(header, VCFHeaderVersion.VCF4_3);
			try(final CloseableIterator<RLine> iter=shuffled.iterator()) {
				while(iter.hasNext())
					{
					final VariantContext ctx= vcfCodec.decode(iter.next().line);
					out.add(ctx);
					}
				}
			shuffled.cleanup();
			shuffled=null;
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally {
			if(shuffled!=null) shuffled.cleanup();
			CloserUtil.close(shuffled);
			}
		}
	
	@Override
	protected int beforeVcf() {
		if(seed==-1L) seed= System.currentTimeMillis();
		return 0;
		}
	

	public static void main(final String[] args)
		{
		new VCFShuffle().instanceMainWithExit(args);
		}
	}
