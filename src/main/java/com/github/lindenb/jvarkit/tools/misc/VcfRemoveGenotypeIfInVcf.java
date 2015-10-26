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

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.vcf.TabixVcfFileReader;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


/**
 * VcfRemoveGenotypeIfInVcf
 * @author lindenb
 *
 */
public class VcfRemoveGenotypeIfInVcf extends AbstractVcfRemoveGenotypeIfInVcf {
	
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(VcfRemoveGenotypeIfInVcf.class);

	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractVcfRemoveGenotypeIfInVcf.AbstractVcfRemoveGenotypeIfInVcfCommand
		{		

	private TabixVcfFileReader tabix=null;
	
	
	@Override
		protected Collection<Throwable> doVcfToVcf(String inputName,
				VcfIterator in, VariantContextWriter out) throws IOException {
		long genotypes_reset_count=0;
		long genotypes_empty=0;
		long variant_changed=0;
		long variant_unchanged=0;
		ArrayList<VariantContext> overlappingList =new ArrayList<VariantContext>();
		VCFHeader h2=in.getHeader();
		addMetaData(h2);
		out.writeHeader(h2);
		while(in.hasNext())
			{
			VariantContext ctx= in.next();
			
			/* no genotypes */
			if(ctx.getNSamples()==0)
				{
				out.add(ctx);
				++variant_unchanged;
				continue;
				}
			/* get overlapping variants */
			overlappingList.clear();
			Iterator<VariantContext> iter = this.tabix.iterator(
					ctx.getContig(), ctx.getStart(), ctx.getEnd()+1);
			while(iter.hasNext())
				{
				VariantContext ctx2=iter.next();
				//check position
				if(ctx2.getStart()!=ctx.getStart() ) continue;
				if(!ctx2.getContig().equals(ctx.getContig()) ) continue;
				//not ALT alleles
				if(ctx2.getAlternateAlleles().isEmpty() ) continue;
				//not same REF anyway
				if(!ctx.getReference().equals(ctx2.getReference())) continue;
				overlappingList.add(ctx2);
				}
			
			/* no overlapping variants */
			if(overlappingList.isEmpty())
				{
				out.add(ctx);
				++variant_unchanged;
				continue;
				}
			
			ArrayList<Genotype> genotypes=new ArrayList<Genotype>(ctx.getGenotypes());
			if(genotypes.size()!=h2.getNGenotypeSamples()) throw new RuntimeException();
		
			int genotype_index=0;
			while(genotype_index< genotypes.size())
				{
				Genotype g=genotypes.get(genotype_index);
				if(g.isNoCall() || !g.isCalled() || g.isHomRef())
					{
					++genotype_index;
					continue;
					}
				//all the alleles for this genotype
				Set<Allele> alt_alleles = new HashSet<Allele>(g.getAlleles());
				//remove REF allele
				alt_alleles.remove(ctx.getReference());
				if(alt_alleles.isEmpty()) throw new IllegalStateException("??");
				
				//remove all the alt alleles found in overlapping
				for(VariantContext ctx2: overlappingList)
					{
					alt_alleles.removeAll(ctx2.getAlternateAlleles());
					}

				
				if(!alt_alleles.isEmpty())/* contains at least one unknown allele */
					{
					++genotype_index;
					}
				else
					{
					genotypes_reset_count++;
					genotypes.remove(genotype_index);
					}
				}
			
			
			//nothing to print
			if(genotypes.isEmpty())
				{
				++genotypes_empty;
				if(removeVariantNoGenotype)
					{
					++variant_changed;
					continue;
					}
				}
			//nothing removed
			if(genotypes.size() == h2.getNGenotypeSamples() )
				{
				out.add(ctx);
				++variant_unchanged;
				continue;
				}
			VariantContextBuilder vcb=new VariantContextBuilder(ctx);
			vcb.genotypes(genotypes);
			out.add(vcb.make());
			++variant_changed;
			}
		LOG.info("Number of genotype set to ./. :"+genotypes_reset_count);
		LOG.info("Number of variants without any ALT genotype:"+genotypes_empty);
		LOG.info("Number of variants changed:"+variant_changed);
		LOG.info("Number of variants unchanged:"+variant_unchanged);
		return RETURN_OK;
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception
			{
			if( this.tabixFilePath==null)
				{
				return wrapException("Undefined VCF tabix  file");
				}
			
			try
				{
				LOG.info("Opening "+tabixFilePath);
				this.tabix=new TabixVcfFileReader(tabixFilePath);
				return doVcfToVcf(inputName);
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(this.tabix);
				this.tabix=null;
				}
			}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new VcfRemoveGenotypeIfInVcf().instanceMainWithExit(args);
	}

}
