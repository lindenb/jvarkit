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
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.PrintStream;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class Biostar130456 extends AbstractBiostar130456
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Biostar130456.class);

	private final static String SAMPLE_TAG="__SAMPLE__";
	
	
		@Override
		public Collection<Throwable> call(final String inputName) throws Exception
			{
			if(super.filepattern==null || !filepattern.contains(SAMPLE_TAG))
				{
				return wrapException("File pattern is missing "+SAMPLE_TAG);
				}
			PrintStream out = null;
			VcfIterator in=null;
			try
				{
				out = openFileOrStdoutAsPrintStream();
				in = super.openVcfIterator(inputName);
				final VCFHeader header=in.getHeader();
				final Set<String> samples = new HashSet<String>(header.getSampleNamesInOrder());
				final Map<String,VariantContextWriter> sample2writer=new HashMap<String,VariantContextWriter>(samples.size());
	
				if(samples.isEmpty())
					{
					return wrapException("VCF doesn't contain any sample");
					}
				LOG.info("N sample:"+samples.size());
				for(final String sample:samples)
					{
					final VCFHeader h2=new VCFHeader(
							header.getMetaDataInInputOrder(),
							Collections.singleton(sample)
							);
					super.addMetaData(h2);
					final String sampleFile= filepattern.replaceAll(SAMPLE_TAG,sample);
					out.println(sampleFile);
					final File fout = new File(sampleFile);
					if(fout.getParentFile()!=null) fout.getParentFile().mkdirs();
					final VariantContextWriter w= VCFUtils.createVariantContextWriter(fout);
					w.writeHeader(h2);
					sample2writer.put(sample, w);
					}
				final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
				while(in.hasNext())
					{
					final VariantContext ctx= progress.watch(in.next());
					for(final String sample: samples)
						{
						final Genotype g= ctx.getGenotype(sample);
						if(g==null) continue;
						if(remove_uncalled && (!g.isAvailable() || !g.isCalled() || g.isNoCall()))
							{
							continue;
							}
						if(remove_homref && g.isHomRef()) continue;
						final VariantContextWriter w= sample2writer.get(sample);
						final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
						final GenotypeBuilder gb=new GenotypeBuilder(g);
						vcb.genotypes(Collections.singletonList(gb.make()));
						final VariantContext ctx2= vcb.make();
						w.add(ctx2);
						}
					}
				for(final String sample:samples)
					{
					LOG.info("Closing for sample "+sample);
					final VariantContextWriter w= sample2writer.get(sample);
					w.close();
					}
				progress.finish();
				out.flush();
				return RETURN_OK;
				}
			catch (Exception e)
				{
				return wrapException(e);
				}
			finally
				{
				CloserUtil.close(out);
				CloserUtil.close(in);
				}
			}
		

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new Biostar130456().instanceMainWithExit(args);
		}

	}
