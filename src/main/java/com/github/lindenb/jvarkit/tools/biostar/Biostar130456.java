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
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class Biostar130456 extends AbstractBiostar130456
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Biostar130456.class);

	private final static String SAMPLE_TAG="__SAMPLE__";
	
	@Override
	public Command createCommand() {
		return new MyCommand();
		}
	
	static private class MyCommand extends AbstractBiostar130456.AbstractBiostar130456Command
		{
		@Override
		public Collection<Throwable> call() throws Exception
			{
			final List<String> args =  getInputFiles();
			
			if(super.filepattern==null || !filepattern.contains(SAMPLE_TAG))
				{
				return wrapException("File pattern is missing "+SAMPLE_TAG);
				}
			
			VcfIterator in=null;
			try
				{
				if(args.isEmpty())
					{
					LOG.info("Reading from <stdin>");
					in = VCFUtils.createVcfIteratorFromInputStream(stdin());
					}
				else if(args.size()==1)
					{
					String filename=args.get(0);
					LOG.info("Reading from "+filename);
					in = VCFUtils.createVcfIterator(filename);
					}
				else
					{
					return wrapException(getMessageBundle("illegal.number.of.arguments"));
					}
				VCFHeader header=in.getHeader();
				Set<String> samples = new HashSet<String>(header.getSampleNamesInOrder());
				Map<String,VariantContextWriter> sample2writer=new HashMap<String,VariantContextWriter>(samples.size());
	
				if(samples.isEmpty())
					{
					return wrapException("VCF doesn't contain any sample");
					}
				LOG.info("N sample:"+samples.size());
				for(String sample:samples)
					{
					
					VCFHeader h2=new VCFHeader(
							header.getMetaDataInInputOrder(),
							Collections.singleton(sample)
							);
					h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
					h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
					h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
					h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
					String sampleFile= filepattern.replaceAll(SAMPLE_TAG,sample);
					stdout().println(sampleFile);
					File fout = new File(sampleFile);
					if(fout.getParentFile()!=null) fout.getParentFile().mkdirs();
					VariantContextWriter w= VCFUtils.createVariantContextWriter(fout);
					w.writeHeader(h2);
					sample2writer.put(sample, w);
					}
				SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
				while(in.hasNext())
					{
					VariantContext ctx= progress.watch(in.next());
					for(String sample: samples)
						{
						Genotype g= ctx.getGenotype(sample);
						if(g==null) continue;
						if(remove_uncalled && (!g.isAvailable() || !g.isCalled() || g.isNoCall()))
							{
							continue;
							}
						if(remove_homref && g.isHomRef()) continue;
						VariantContextWriter w= sample2writer.get(sample);
						VariantContextBuilder vcb=new VariantContextBuilder(ctx);
						GenotypeBuilder gb=new GenotypeBuilder(g);
						vcb.genotypes(Collections.singletonList(gb.make()));
						VariantContext ctx2= vcb.make();
						w.add(ctx2);
						}
					}
				for(String sample:samples)
					{
					LOG.info("Closing for sample "+sample);
					VariantContextWriter w= sample2writer.get(sample);
					w.close();
					}
				progress.finish();
				return Collections.emptyList();
				}
			catch (Exception e)
				{
				return wrapException(e);
				}
			finally
				{
				CloserUtil.close(in);
				}
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
