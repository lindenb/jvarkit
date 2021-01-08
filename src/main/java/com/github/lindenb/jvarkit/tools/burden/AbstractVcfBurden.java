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

*/
package com.github.lindenb.jvarkit.tools.burden;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.pedigree.Sample;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFReader;


public abstract class AbstractVcfBurden
extends Launcher
{
	private static final Logger LOG = Logger.build(AbstractVcfBurden.class).make();
	protected static final String BURDEN_KEY = "BURDEN_KEY";
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-p","--ped","--pedigree"},description=PedigreeParser.OPT_DESC,required=true)
	private Path pedFile = null;
	@Parameter(names={"-save-vcf","--save-vcf"},description="Save Matching variants into that VCF.")
	private Path outputVcfPath = null;

	
	private Pedigree pedigree = null;
	private Set<Sample> cases = null;
	private Set<Sample> controls = null;
	
	private enum SuperVariant
		{
		SV0,AT_LEAST_ONE_VARIANT
		}
	
	
	protected static class FisherResult {
		double p_value = 0;
		int affected_alt = 0;
		int affected_hom = 0;
		int unaffected_alt = 0;
		int unaffected_hom = 0;

	}
	
	protected FisherResult runFisher(
		final List<VariantContext> variants	
		) {
		final Map<Sample,SuperVariant> indi2supervariant = new HashMap<>(this.cases.size() + this.controls.size());
		this.cases.stream().forEach(S->indi2supervariant.put(S,SuperVariant.SV0));
		this.controls.stream().forEach(S->indi2supervariant.put(S,SuperVariant.SV0));
		
		for(final VariantContext ctx:variants) {
			final Allele observed_alt = ctx.getAltAlleleWithHighestAlleleCount();
			for(final Sample sample : indi2supervariant.keySet() ) {
				if(indi2supervariant.get(sample)==SuperVariant.AT_LEAST_ONE_VARIANT) continue;
				final Genotype g = ctx.getGenotype(sample.getId());	
				if(g==null || g.isHomRef() || g.isNoCall()) continue;
				if( g.getAlleles().stream().anyMatch(A->A.equals(observed_alt))) {
					indi2supervariant.put(sample,SuperVariant.AT_LEAST_ONE_VARIANT);
					break;
					}
				}// end for sample
			}//end of forVariant
		
		final FisherResult fisher = new FisherResult();
		
		for(final Sample sample : indi2supervariant.keySet() ) {
			final SuperVariant superVariant = indi2supervariant.get(sample);
			if(superVariant==SuperVariant.SV0 ) {
				if(sample.isAffected()) fisher.affected_hom++;
				else fisher.unaffected_hom++;
				}
			else // AT_LEAST_ONE_VARIANT 
				{
				if(sample.isAffected()) fisher.affected_alt++;
				else fisher.unaffected_alt++;
				}
			}//end of sample

		
		final FisherExactTest test = FisherExactTest.compute(
				fisher.affected_alt, fisher.affected_hom, 
				fisher.unaffected_alt, fisher.unaffected_hom
				);
		fisher.p_value = test.getAsDouble();
		return fisher;
		}

	
	protected abstract void runBurden(
			final PrintWriter out,
			final VCFReader vcfReader,
			final VariantContextWriter vcw
			) throws IOException;
	
	@Override
	public int doWork(final List<String> args) {
		PrintWriter pw = null;
		VCFReader vcfReader = null;
		VariantContextWriter vcw = null;
		try {
			this.pedigree = new PedigreeParser().parse(this.pedFile);
			
			this.cases = new HashSet<>(this.pedigree.getAffectedSamples());
			this.controls = new HashSet<>(this.pedigree.getUnaffectedSamples());
			
			
			final String vcfIn = super.oneAndOnlyOneFile(args);
			vcfReader = VCFReaderFactory.makeDefault().open(Paths.get(vcfIn),true);
			final VCFHeader header = vcfReader.getHeader();
			final Set<String> samplesInVcf = new HashSet<>(header.getSampleNamesInOrder());
			
			if(this.outputVcfPath!=null) {
				vcw = VCFUtils.createVariantContextWriterToPath(this.outputVcfPath);
				header.addMetaDataLine(new VCFInfoHeaderLine(BURDEN_KEY, 1, VCFHeaderLineType.String,"Burden key"));
				JVarkitVersion.getInstance().addMetaData(this, header);
				vcw.writeHeader(header);
			}
			
			
			this.cases.removeIf(S->!samplesInVcf.contains(S.getId()));
			this.controls.removeIf(S->!samplesInVcf.contains(S.getId()));
			
			if(this.cases.isEmpty()) {
				LOG.error("no affected in "+this.pedFile);
				return -1;
				}
			if(this.controls.isEmpty()) {
				LOG.error("no controls in "+this.pedFile);
				return -1;
				}
			pw = super.openPathOrStdoutAsPrintWriter(this.outputFile);
			
			runBurden(pw,vcfReader,vcw);
			
			pw.flush();
			pw.close();
			pw=null;
			vcfReader.close();
			vcfReader = null;
			return 0;
			} 
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally {
			CloserUtil.close(vcfReader);
			CloserUtil.close(pw);
			CloserUtil.close(vcw);
			}
	
		}
	
}
