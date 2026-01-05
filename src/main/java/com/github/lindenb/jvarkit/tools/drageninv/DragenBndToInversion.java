/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.drageninv;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.dict.OrderChecker;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.AutoMap;
import com.github.lindenb.jvarkit.util.JVarkitVersion;

import htsjdk.samtools.util.CoordMath;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

# Motivation 

Dragen has no VCF with SVTYPE=INV

[https://jp.support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/MantaInversions_fDG.htm](https://jp.support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/MantaInversions_fDG.htm)

> Inversions are reported as a set of breakends. For example, given a simple reciprocal inversion, four breakends are reported, sharing the same EVENT INFO tag. The following is an example breakend records representing a simple reciprocal inversion:

END_DOC
 */
@Program(
		name="dragenbnd2inv",
		description="Converts Dragen BND to inversions",
		keywords={"vcf","dragen","inversion","bnd","sv"},
		creationDate="20241016",
		modificationDate="20241016")
public class DragenBndToInversion extends OnePassVcfLauncher {
	private static final String EVENT_KEY ="EVENT";
	private static Logger LOG=Logger.of(DragenBndToInversion.class);
	@Parameter(names={"--emit-other"},description="write other variants that are not identified as BND/INV")
	private boolean emit_other = false;
	@Parameter(names={"--emit-original"},description="emit original BND variants that makes an NV")
	private boolean emit_original = false;
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iter, VariantContextWriter out) {
		final VCFHeader header= iter.getHeader();
		if(header.getInfoHeaderLine(EVENT_KEY)==null) {
			LOG.error("INFO/"+EVENT_KEY+" is not defined in "+inputName);
			return -1;
			}
		if(header.getInfoHeaderLine("SVLEN")==null) {
			header.addMetaDataLine(new VCFInfoHeaderLine("SVLEN", 1, VCFHeaderLineType.Integer,"SV length"));
			return -1;
			}
		if(header.getNGenotypeSamples()>1) {
			LOG.error("Cannot handle more than one genotyped sample");
			return -1;
			}
		JVarkitVersion.getInstance().addMetaData(this, header);
		out.writeHeader(header);
		final AutoMap<String, VariantContext,List<VariantContext>> event2variants = AutoMap.makeList();
		String prev_chrom="";
		final OrderChecker<VariantContext> checker= new OrderChecker<>();
		final List<Allele> alleles= Arrays.asList(Allele.REF_N,Allele.create("<INV>", false));

		for(;;) {
			final VariantContext ctx0 = iter.hasNext()? checker.apply(iter.next()) : null;
			if(ctx0==null || !ctx0.getContig().equals(prev_chrom)) {
				for(String eventid:event2variants.keySet()) {
					final List<VariantContext> variants = event2variants.get(eventid).
							stream().
							sorted((A,B)->Integer.compare(A.getStart(),B.getStart())).
							collect(Collectors.toList());
					// should be on the same contig
					if(variants.size()!=4 || variants.stream().map(V->V.getContig()).collect(Collectors.toSet()).size()!=1) {
						if(this.emit_other) {
							for(VariantContext ctx2:variants) {
								out.add(ctx2);
								}
							}
						continue;
						}
					// different genotypes ??
					if(header.getNGenotypeSamples()==1 && variants.stream().map(V->V.getGenotype(0).getType()).collect(Collectors.toSet()).size()!=1) {
						LOG.warn("different genotypes for INFO/"+EVENT_KEY+"="+eventid);
						if(this.emit_other) {
							for(VariantContext ctx2:variants) {
								out.add(ctx2);
								}
							}
						continue;
						}
					final VariantContextBuilder vcb=new VariantContextBuilder(
						inputName,
						variants.get(0).getContig(),
						variants.get(0).getStart(),
						variants.get(3).getEnd(),
						alleles
						);
					vcb.id(variants.get(0).getID());
					vcb.log10PError(variants.get(0).getLog10PError());
					vcb.attribute(VCFConstants.END_KEY,variants.get(3).getEnd());
					vcb.attribute("SVLEN",CoordMath.getLength(variants.get(0).getStart(), variants.get(3).getEnd()));
					vcb.attribute(VCFConstants.SVTYPE,"INV");
					vcb.attribute(EVENT_KEY,eventid);
					/*
					if(variants.stream().allMatch(V->V.hasAttribute("CIPOS"))) {
						int min_ci5 = variants.stream().mapToInt(V->V.getStart()+V.getAttributeAsIntList("CIPOS", 0).get(0)).max().orElse(variants.get(0).getStart()) - variants.get(0).getStart();
						int max_ci3 = variants.get(3).getStart() - variants.stream().mapToInt(V->V.getStart()+V.getAttributeAsIntList("CIPOS", 0).get(1)).max().orElse(variants.get(0).getStart());
						vcb.attribute("CIPOS", Arrays.asList(min_ci5,max_ci3));
						}*/
					if(variants.get(0).hasAttribute("JUNCTION_QUAL")) {
						vcb.attribute("JUNCTION_QUAL",variants.get(0).getAttribute("JUNCTION_QUAL"));
						}
					if(header.getNGenotypeSamples()==1) {
						final GenotypeBuilder gb=new GenotypeBuilder(variants.get(0).getGenotype(0));
						gb.alleles(variants.get(0).getGenotype(0).getAlleles().stream().map(A->A.isReference()?alleles.get(0):alleles.get(1)).collect(Collectors.toList()));
						vcb.genotypes(Arrays.asList(gb.make()));
						}
					
					out.add(vcb.make());
					if(emit_original) {
						for(VariantContext ctx2:variants) {
							out.add(ctx2);
							}
						}
					}
				if(ctx0==null)  break;
				prev_chrom = ctx0.getContig();
				event2variants.clear();
				}
			final String svType = ctx0.getAttributeAsString(VCFConstants.SVTYPE,"");
			if(!svType.equals("BND") || !ctx0.hasAttribute(EVENT_KEY) ) {
				if(this.emit_other) out.add(ctx0);
				continue;
				}
			final String event = ctx0.getAttributeAsString(EVENT_KEY, "");
			event2variants.insert(event, ctx0);
			}
		return 0;
		}
	
	
public static void main(final String[] args) {
	new DragenBndToInversion().instanceMainWithExit(args);
	}
}
