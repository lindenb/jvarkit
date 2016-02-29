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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import javax.script.Bindings;
import javax.script.CompiledScript;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


/**
 * VcfRemoveGenotypeJs
 * @author lindenb
 *
 */
public class VcfRemoveGenotypeJs extends AbstractVcfRemoveGenotypeJs {
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfRemoveGenotypeJs.class);

	private CompiledScript  script=null;

	@Override
	protected Collection<Throwable> doVcfToVcf(String inputName, VcfIterator in, VariantContextWriter out)
			throws IOException {
		try {
		this.script  = super.compileJavascript();
		final	VCFHeader h2=new VCFHeader(in.getHeader());
		addMetaData(h2);
		out.writeHeader(h2);
		final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(in.getHeader());
        final Bindings bindings = this.script.getEngine().createBindings();
        bindings.put("header", in.getHeader());

		while(in.hasNext())
			{
			final VariantContext ctx = in.next();
			bindings.put("variant", ctx);
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			List<Genotype> genotypes= new ArrayList<>();
			int countCalled = ctx.getNSamples();
			
			for(int i=0;i< ctx.getNSamples();++i)
				{
				Genotype genotype = ctx.getGenotype(i);
				bindings.put("genotype", genotype);
				
				if(genotype.isNoCall() || !genotype.isAvailable())
					{
					countCalled--;
					}
				else if(genotype.isCalled() &&
					!super.evalJavaScriptBoolean(this.script, bindings)) {
					if(super.replaceByHomRef){
						List<Allele> homRefList=new ArrayList<>(genotype.getPloidy());
						for(int p=0;p< genotype.getPloidy();++p)
							{
							homRefList.add(ctx.getReference());
							}
						genotype = new GenotypeBuilder(genotype).alleles(homRefList).make();
						} 
					else
						{
						genotype = GenotypeBuilder.createMissing(genotype.getSampleName(), genotype.getPloidy());
						}
					countCalled--;
					}
				genotypes.add(genotype);
				}
			if(countCalled==0 && super.removeCtxNoGenotype) {
				continue;
			}
			vcb.genotypes(genotypes);
			out.add(vcb.make());
			}
		
		progress.finish();
		
		return RETURN_OK;
		} catch(Exception err) {
			return wrapException(err);
		} finally
		{
			this.script=null;
		}
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		return doVcfToVcf(inputName);
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new VcfRemoveGenotypeJs().instanceMainWithExit(args);
	}

}
