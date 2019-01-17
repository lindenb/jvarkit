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
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import com.github.lindenb.jvarkit.util.vcf.TabixVcfFileReader;
import htsjdk.variant.vcf.VCFIterator;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;


/**
BEGIN_DOC
## Example

```bash
$ gunzip -c input.vcf.gz | grep -E '1308871'
1	1308871	.	A	T	6.2	.	.	GT:PL:DP:GQ	1/1:35,3,0:1:4



$ gunzip -c input.vcf.gz | grep -E '(^#|1308871)' |\
  java -jar dist/vcfresetvcf.jar -x ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz |\
  grep -v '^#'
1	1308871	.	A	T	6.20	.	.	GT	./.

```
END_DOC

 */
@Program(name="vcfresetvcf",
	description="Reset Genotypes in VCF (./.) if they've been found in another VCF indexed with tabix",
	keywords={"vcf","genotype"}
)
public class VcfRemoveGenotypeIfInVcf extends Launcher {
	private TabixVcfFileReader tabix=null;
	
	@Parameter(names="-x",description="remove variant if there is no called genotype")
	private boolean removeVariantNoGenotype=false;
	

	private static final Logger LOG = Logger.build(VcfRemoveGenotypeIfInVcf.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator in, VariantContextWriter out) {
		long genotypes_reset_count=0;
		long genotypes_empty=0;
		long variant_changed=0;
		long variant_unchanged=0;
		final ArrayList<VariantContext> overlappingList =new ArrayList<VariantContext>();
		final VCFHeader h2=in.getHeader();
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		JVarkitVersion.getInstance().addMetaData(getClass().getSimpleName(), h2);
		out.writeHeader(h2);
		while(in.hasNext())
			{
			final VariantContext ctx= in.next();
			
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
		return 0;
		}
	
	@Parameter(names={"-t","--tabix"},description="Tabix indexed VCF file",required=true)
	private String tabixFilePath=null;
	
	@Override
	public int doWork(List<String> args) {

		if( tabixFilePath==null)
			{
			LOG.error("Undefined VCF tabix  file");
			return -1;
			}
		
		try
			{
			LOG.info("Opening "+tabixFilePath);
			this.tabix=new TabixVcfFileReader(tabixFilePath);
			return doVcfToVcf(args, outputFile);
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.tabix);
			}
		}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new VcfRemoveGenotypeIfInVcf().instanceMainWithExit(args);
	}

}
