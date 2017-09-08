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
* 2015 moved to htsjdk + knime

*/
package com.github.lindenb.jvarkit.tools.vcfbigwig;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
/*
BEGIN_DOC

## Example

```bash
 java -jar dist/vcfbigwig.jar \
 	-T GERP \
 	-B gerp.bw input.vcf.gz 
	
##INFO=<ID=GERP,Number=1,Type=Float,Description="Values from bigwig file: com.github.lindenb.jvarkit.tools.vcfbigwig.VCFBigWig BIGWIG=gerp.bw TAG=GERP IN=input.vcf.gz    VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO(...)
A	33926	.	G	A	182	.	GERP=-6.35(...)
A	45365	.	A	G	222	.	GERP=-3.55(...)
```


END_DOC
*/
@Program(name="vcfbigwig",
	description="Annotate a VCF with values from a bigwig file",
	keywords={"vcf","wig","wiggle","bigwig"}
	)
public class VCFBigWig extends Launcher
	{

	private static final Logger LOG = Logger.build(VCFBigWig.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();

	
	public static class CtxWriterFactory 
	implements VariantContextWriterFactory
		{
		@Parameter(names={"-B","--bigwig"},description="Path to the bigwig file",required=true)
		private String biwWigFile = null;
	
		@Parameter(names={"-T","--tag","-tag"},description="Name of the INFO tag. default: name of the bigwig")
		private String TAG = null;
	
		@Parameter(names={"-C","--contained"},description="Specifies wig values must be contained by region. if false: return any intersecting region values")
		private boolean contained = false;
	
		@Parameter(names={"-a","--aggregate"},description="How to aggregate overlapping values: 'avg' average; 'median': median, 'first': use first, 'all' : print all the data")
		private AggregateMethod aggregateMethod  = AggregateMethod.avg;
	
		@Parameter(names={"-t","--transform"},description="Deprecated",hidden=true)
		private String _convertChrName = null;
	
		@Parameter(names={"--onNotFound"},description="[20170707] " + ContigNameConverter.OPT_ON_NT_FOUND_DESC)
		private ContigNameConverter.OnNotFound onContigNotFound =ContigNameConverter.OnNotFound.SKIP;
		
		
		private BBFileReader bbFileReader=null;

		
		private class CtxWriter extends DelegateVariantContextWriter
			{
			private final AggregateMethod aggregateMethod;
			private final ContigNameConverter contigNameConverter;
			private final Set<String> userContigsNotFound = new HashSet<>();
			private final List<Float> values=new ArrayList<Float>();

			
			CtxWriter(final VariantContextWriter delegate) {
				super(delegate);
				this.aggregateMethod = CtxWriterFactory.this.aggregateMethod;
				this.contigNameConverter = ContigNameConverter.fromContigSet(new HashSet<>(CtxWriterFactory.this.bbFileReader.getChromosomeNames()));
				this.contigNameConverter.setOnNotFound(CtxWriterFactory.this.onContigNotFound);
				}
			
			
			@Override
			public void writeHeader(final VCFHeader header) {
					
				final VCFHeader h2=new VCFHeader(header);
				
				if(h2.getInfoHeaderLine(CtxWriterFactory.this.TAG)!=null)
					{
					throw new JvarkitException.DuplicateVcfHeaderInfo(h2,CtxWriterFactory.this.TAG);
					}
				
				if(this.aggregateMethod.equals(AggregateMethod.all))
					{
					h2.addMetaDataLine(new VCFInfoHeaderLine(
							CtxWriterFactory.this.TAG,
							VCFHeaderLineCount.UNBOUNDED,
							VCFHeaderLineType.Float,
							"Values from bigwig file: "+CtxWriterFactory.this.biwWigFile
							));
					}
				else
					{
					h2.addMetaDataLine(new VCFInfoHeaderLine(
							CtxWriterFactory.this.TAG,1,
							VCFHeaderLineType.Float,
							"Values from bigwig file: "+CtxWriterFactory.this.biwWigFile
							));
					}
				
				super.writeHeader(h2);
				}

			@Override
			public void add(final VariantContext ctx) {
				this.values.clear();
				final String variantChrom=  contigNameConverter.apply(ctx.getContig());
				
				if( variantChrom == null) {
					if(!this.userContigsNotFound.contains(ctx.getContig()))
						{
						this.userContigsNotFound.add(ctx.getContig());
						LOG.warn("Bigwig file \""+CtxWriterFactory.this.biwWigFile+"\" doesn't contains contig "+ variantChrom+"/"+ctx.getContig());
						}
					super.add(ctx);
					return;
					}
				
				
				final BigWigIterator iter= CtxWriterFactory.this.bbFileReader.getBigWigIterator(
						variantChrom,
						ctx.getStart()-1,
						variantChrom,
						ctx.getEnd(),
						CtxWriterFactory.this.contained
						);
				while(iter!=null && iter.hasNext())
					{
					final WigItem item=iter.next();
					final float v=item.getWigValue();
					values.add(v);
					if(this.aggregateMethod.equals(AggregateMethod.first)) break;
					}
				
				if(values.isEmpty())
					{
					super.add(ctx);
					return;
					}
				final VariantContextBuilder b=new VariantContextBuilder(ctx);

				switch(this.aggregateMethod)
					{
					case all:
						b.attribute(CtxWriterFactory.this.TAG,values);
						break;
					case avg:
						double total=0L;
						for(final Float f:values) total+=f;
						b.attribute(CtxWriterFactory.this.TAG,(float)(total/values.size()));
						break;
					case first:
						b.attribute(CtxWriterFactory.this.TAG,values.get(0));
						break;
					case median:
						final double median_value;
						values.sort((A,B)->A.compareTo(B));
						final int mid_x= values.size()/2;
						if(values.size()==1)
							{
							median_value = values.get(0);
							}
						else if(values.size()%2==0)
	                        {
	                		median_value =  (values.get(mid_x-1)+values.get(mid_x))/2.0;
	                        }
		                else
	                        {
	                		median_value =  values.get(mid_x);
	                        }
						b.attribute(CtxWriterFactory.this.TAG,median_value);
						break;
					default: throw new IllegalStateException();
					}
				super.add(b.make());
				}
			
			@Override
			public void close() {
				if(!this.userContigsNotFound.isEmpty())
					{
					LOG.warn("\""+CtxWriterFactory.this.biwWigFile+
							"\": Contigs not found :"+
							String.join(" ", userContigsNotFound));
					}
				super.close();
				}
			
			}
		
		@Override
		public int initialize() {
			if(this.biwWigFile==null || this.biwWigFile.isEmpty())
				{
				LOG.info("Undefined BigWig file ");
				return -1;
				}
		
			try
				{
				this.bbFileReader= new BBFileReader(this.biwWigFile);
				if(!this.bbFileReader.isBigWigFile())
					{
					this.bbFileReader=null;
					throw new IOException(this.biwWigFile+" is not a bigWIG file.");
					}
		
				if(this.TAG==null || this.TAG.isEmpty())
					{
					this.TAG=this.biwWigFile;
					int i=TAG.lastIndexOf(File.separator);
					if(i!=-1) TAG=TAG.substring(i+1);
					i=this.TAG.indexOf('.');
					this.TAG=this.TAG.substring(0,i);
					LOG.info("setting tag to "+this.TAG);
					}
				return 0;
				}
			catch(final Exception err)
				{
				LOG.error(err);
				return -1;
				}
			}
		
		@Override
		public VariantContextWriter open(VariantContextWriter delegate) {
			return new CtxWriter(delegate);
			}

		@Override
		public void close() throws IOException {
			try
				{
				if(this.bbFileReader!=null)
					{
					CloserUtil.close(this.bbFileReader.getBBFis());
					}
				CloserUtil.close(this.bbFileReader);
				this.bbFileReader=null;
				}
			catch(final Exception err)
				{
				LOG.error("Error",err);
				}
			VariantContextWriterFactory.super.close();
			}
		}
	
	
	private enum AggregateMethod
		{
		avg,median,first,all
		}
	public VCFBigWig()
		{
		}
	
	
	@Override
	protected int doVcfToVcf(final String inputName, final VcfIterator r, final VariantContextWriter delegate) {
		final VariantContextWriter w = this.component.open(delegate);
		w.writeHeader(r.getHeader());
		final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(r.getHeader()).logger(LOG);
		while(r.hasNext())
			{
			w.add(progress.watch(r.next()));
			
			// JVM crash sometimes ? suspect there is a memory leak ?
			if(progress.getCount()%1000L==0)
				{
				System.gc();
				}
			}
		progress.finish();
		w.close();
		return 0;
		}
	
	@Override
	public int doWork(final List<String> args) {
		try 
			{
			if(this.component.initialize()!=0) {
				return -1;
				}
			return doVcfToVcf(args,this.outputFile);
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.component);
			}
		}
	
	public static void main(final String[] args) throws IOException
		{
		new VCFBigWig().instanceMainWithExit(args);
		}
}
