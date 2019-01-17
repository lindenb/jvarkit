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
package com.github.lindenb.jvarkit.tools.skat;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlType;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VCFBuffer;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
/**
BEGIN_DOC

## Example

```
$ java -jar dist/vcfskat.jar vcf_with_samples.vcf | grep SKAT
##SKAT=0.2215079
```

```
$ java -jar dist/vcfskat.jar -p vcf_with_samples.vcf 
0.2215079
```


END_DOC

 */
@Program(
		name="vcfskat",
		description="Calculate SKAT score for a VCF.",
		keywords={"vcf","pedigree","skat","burden"}
		)
public class VcfSkat extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfSkat.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-p","--print"},description="Don't print the VCF on output,Just print the score on ouput, don't print the VCF itself")
	private boolean just_print_skat_result=false;
	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();

	
	@XmlType(name="vcfskat")
	@XmlRootElement(name="vcfskat")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class CtxWriterFactory 
	implements VariantContextWriterFactory
		{
		@ParametersDelegate
		private SkatFactory skat = new SkatFactory();
		@Parameter(names={"-ped","--pedigree"},description=Pedigree.OPT_DESCRIPTION+" If not defined, I will try to extract the pedigree from  the VCFheader.")
		private File pedigreeFile=null;
		@Parameter(names={"-tag","--tag"},description="tag name in VCF header \"tag=p-value\"")
		private String tagName = "SKAT";
		
		private class CtxWriter extends DelegateVariantContextWriter
			{
			private final SkatFactory.SkatExecutor executor = CtxWriterFactory.this.skat.build();
			private final File pedigreeFile = CtxWriterFactory.this.pedigreeFile;
			private final String tagName = CtxWriterFactory.this.tagName;
			private Double p_value = null;
			
			final VCFBuffer buffer= new VCFBuffer();
			CtxWriter(final VariantContextWriter delegate) {
				super(delegate);
				}
			@Override
			public void writeHeader(final VCFHeader header) {
				this.buffer.writeHeader(header);
				}
			@Override
			public void add(final VariantContext ctx) {
				buffer.add(ctx);
				}
			@Override
			public void close() {
				this.buffer.close();
				
				final Pedigree pedigree;
				if(this.pedigreeFile!=null)
					{
					try {
						pedigree = new Pedigree.Parser().parse(this.pedigreeFile);
						}
					catch(final IOException err)
						{
						throw new RuntimeIOException(err);
						}
					}
				else
					{
					pedigree = new Pedigree.Parser().parse(this.buffer.getHeader());
					}
				final Set<Pedigree.Person> samples= new HashSet<>(
						pedigree.getPersons()
						);
				
				
				final List<VariantContext> variants = this.buffer.stream().
					filter(this.executor.getUpstreamVariantFilter()).
					collect(Collectors.toList()
					);
				final VCFHeader h2 = new VCFHeader(buffer.getHeader());
				final SkatFactory.SkatResult result=  this.executor.execute(variants, samples);
				if(result.isError())
					{
					LOG.warn(result.getMessage());
					h2.addMetaDataLine(new VCFHeaderLine(this.tagName+".error",result.getMessage()));
					}
				else
					{
					this.p_value = result.getPValue();
					LOG.info(this.tagName +"="+this.p_value);
					h2.addMetaDataLine(new VCFHeaderLine(this.tagName, String.valueOf(this.p_value)));
					}
				super.writeHeader(h2);
				this.buffer.stream().forEach(V->super.add(V));
				this.buffer.dispose();
				super.close();
				}
			}
		@Override
		public CtxWriter open(VariantContextWriter delegate) {
			return new CtxWriter(delegate);
			}
		}
	
	@Override
	public int doVcfToVcf(final String inputName,final VCFIterator r, final VariantContextWriter delegate)
		{
		final CtxWriterFactory.CtxWriter out = this.component.open(delegate);
		final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(r.getHeader()).logger(LOG);
		out.writeHeader(r.getHeader());
		while(r.hasNext())
			{
			final VariantContext ctx = progress.watch(r.next());
			if(this.just_print_skat_result)
				{
				if(!out.executor.getUpstreamVariantFilter().test(ctx)) continue;
				}
			out.add(ctx);
			}
		out.close();
		progress.finish();
		if(this.just_print_skat_result)
			{
			PrintWriter pw=null;
			try {
				pw = super.openFileOrStdoutAsPrintWriter(this.outputFile);
				if(out.p_value!=null) pw.println(out.p_value);
				pw.flush();
				pw.close();
				pw=null;
				return 0;
			} catch(final IOException err)
				{
				LOG.error(err);
				return -1;
				}
			finally
				{
				CloserUtil.close(pw);
				}
			}
		return 0;
		}
	
	@Override
	protected VariantContextWriter openVariantContextWriter(File outorNull) throws IOException {
	if(this.just_print_skat_result)
		{
		return VCFUtils.createVariantContextWriterToOutputStream(new NullOuputStream());
		}
	return super.openVariantContextWriter(outorNull);
	}

		
	@Override
	public int doWork(final List<String> args)
		{
		try {
			if(this.component.initialize()!=0) return -1;
			return doVcfToVcf(args, this.outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.component);
			}
		}
	
	public static void main(final String[] args) {
		new VcfSkat().instanceMainWithExit(args);
	}
	}
