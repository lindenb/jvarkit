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
package com.github.lindenb.jvarkit.tools.vcfallele2symbolic;



import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/**
BEGIN_DOC

## Motivation

for SV vcf alleles, the annotation with VEP can be heavy if the alleles are coded as loarge IUPAC string

## Example


```bash
java -jar dist/jvarkit.jar vcfallele2symbolic   in.vcf > out.vcf
```

END_DOC
 */
@Program(name="vcfallele2symbolic",
	description="Convert large IUPAC allele to symbolic",
	keywords={"vcf","allele","snv","vcf"},
	creationDate="20250717",
	modificationDate="20250717",
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)
public class VcfAlleleToSymbolic extends OnePassVcfLauncher
	{
	private static final Logger LOG = Logger.of(VcfAlleleToSymbolic.class);

	@Parameter(names={"-n","--length"},description="Change allele if the length it greater or equal to that value")
	private int allele_length = 10;
	@Parameter(names={"--keep-ref"},description="do not replace the REF allele")
	private boolean keep_ref_allele = false;

	public VcfAlleleToSymbolic() {
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	private boolean isBadAllele(final Allele a) {
		return AcidNucleics.isATGCN(a) && a.length()> this.allele_length;
	}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iterin, VariantContextWriter out) {
		final VCFHeader header= iterin.getHeader();
		final Set<VCFHeaderLine> headerLines = new HashSet<>();
		// https://en.wikipedia.org/wiki/Pseudoautosomal_region#Location
		if(header.getInfoHeaderLine(VCFConstants.END_KEY)==null) {
			VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.END_KEY);
			}
		
		
		headerLines.forEach(L->header.addMetaDataLine(L));
		
		JVarkitVersion.getInstance().addMetaData(this, header);
		out.writeHeader(header);
		while(iterin.hasNext()) {
			final VariantContext ctx = iterin.next();
			
			if(ctx.getAlleles().stream().anyMatch(A->isBadAllele(A))) {
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				final String svType = ctx.getAttributeAsString(VCFConstants.SVTYPE, "SYMBOLIC");
				if(!ctx.hasAttribute(VCFConstants.END_KEY)) {
					vcb.attribute(VCFConstants.END_KEY, ctx.getEnd());
					}
				final List<Allele> oldAlleles = ctx.getAlleles();
				final List<Allele> newAlleles = new ArrayList<>(oldAlleles.size());
				for(int i=0;i< oldAlleles.size();i++) {
					Allele a = oldAlleles.get(i);
					if(isBadAllele(a)) {
						if(i==0) {
							// a = Allele.create("<REF>",true); //no , cannot use a symbolic as REF in htsjdk
							if(!keep_ref_allele) {
								a = Allele.REF_N;
								}
							}
						else
							{
							a = Allele.create("<"+svType+(oldAlleles.size()==2?"":String.valueOf(i))+">",false);
							}
						}
					
					newAlleles.add(a);
					}
				vcb.alleles(newAlleles);
				
				if(ctx.hasGenotypes()) {
					final List<Genotype> genotypes= new ArrayList<>(ctx.getNSamples());
					for(Genotype gt: ctx.getGenotypes()) {
						final List<Allele> newGTAlleles = new ArrayList<>(gt.getPloidy());
						for(Allele ga: gt.getAlleles()) {
							if(isBadAllele(ga)) {
								final int i = oldAlleles.indexOf(ga);
								if(i==-1) throw new IllegalArgumentException("cannot find allele "+ga+" in "+oldAlleles);
								newGTAlleles.add(newAlleles.get(i));
								}
							else
								{
								newGTAlleles.add(ga);
								}
							}
						genotypes.add(new GenotypeBuilder(gt).alleles(newGTAlleles).make());
						}
					vcb.genotypes(genotypes);
					}
				out.add(vcb.make());
				}
			else
				{
				out.add(ctx);
				}
			}
		return 0;
		}
	


	public static void main(final String[] args)
		{
		new VcfAlleleToSymbolic().instanceMainWithExit(args);
		}
	
	}
