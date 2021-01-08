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
package com.github.lindenb.jvarkit.tools.vcftrios;

import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BiPredicate;


import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.vcf.JexlGenotypePredicate;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.pedigree.Family;
import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
BEGIN_DOC

## Example


 
END_DOC

 */
@Program(
		name="vcffamilies",
		description="Fills family-based informations in a VCF.",
		keywords={"vcf","pedigree"},
		modificationDate="20200724",
		creationDate="20171123",
		generate_doc=false
		)
public class VCFFamilies
	extends OnePassVcfLauncher
	{
	private static final  Logger LOG = Logger.build(VCFFamilies.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-p","--ped","--pedigree"},description="Pedigree file. "+PedigreeParser.OPT_DESC,required=true)
	private Path pedigreeFile = null;
	@Parameter(names={"-pfx","--prefix"},description="VCF header Attribute prefix")
	private String prefix = "";
	@Parameter(names={"-gf","--genotype-filter"},description=JexlGenotypePredicate.PARAMETER_DESCRIPTION,converter=JexlGenotypePredicate.Converter.class)
	private BiPredicate<VariantContext,Genotype> genotypeFilter = JexlGenotypePredicate.create("!g.isFiltered()");
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();

	

	private class FamilyInfo
		{
		final Family pedFamily;
		private final Set<String> samples=new HashSet<>();
		private final VCFInfoHeaderLine ac;
		private final VCFInfoHeaderLine an;
		private final VCFInfoHeaderLine af;
		
		
		
		FamilyInfo(final Family pedFamily)
			{
			this.pedFamily=pedFamily;
			final String prefix=
					(
					StringUtil.isBlank(VCFFamilies.this.prefix)?
					VCFFamilies.this.prefix+"_":""
					)+
				this.pedFamily.getId()
				;
			this.ac = new VCFInfoHeaderLine(
					prefix+"_" + VCFConstants.ALLELE_COUNT_KEY,
					VCFHeaderLineCount.A,
					VCFHeaderLineType.Integer,
					"Number of alleles carrying the variant allele for the family "+this.pedFamily.getId()
					);
			this.an = new VCFInfoHeaderLine(
					prefix+"_" + VCFConstants.ALLELE_NUMBER_KEY,
					1,
					VCFHeaderLineType.Integer,
					"Number of alleles for the family "+this.pedFamily.getId()
					);
			this.af = new VCFInfoHeaderLine(
					prefix+"_" + VCFConstants.ALLELE_FREQUENCY_KEY,
					VCFHeaderLineCount.A,
					VCFHeaderLineType.Float,
					"Allele Frequency for the family "+this.pedFamily.getId()
					);
			}
		List<VCFInfoHeaderLine> getMetaDataLines() {
			return Arrays.asList(this.ac,this.an,this.af);
			}
		void visit(
				final VariantContextBuilder vcb,
				final VariantContext ctx,
				final List<Allele> alts
				) {
			final int van = 
					this.samples.stream().
					map(S->ctx.getGenotype(S)).
					filter(G->VCFFamilies.this.genotypeFilter.test(ctx, G)).
					mapToInt(G->G.getAlleles().size()).
					sum()
					;
			vcb.attribute(this.an.getID(), van);
			
			
			final int vac[]=new int[alts.size()];
			for(int i=0;i< alts.size();i++)
				{
				final Allele alt = alts.get(i);
				vac[i] = (int)
						this.samples.stream().
						map(S->ctx.getGenotype(S)).
						filter(G->VCFFamilies.this.genotypeFilter.test(ctx, G)).
						flatMap(G->G.getAlleles().stream()).
						filter(A->A.equals(alt)).
						count()
						;
				
				}
			vcb.attribute(this.ac.getID(), vac);
			
			if(van>0)
				{
				final double afs[]=Arrays.stream(vac).mapToDouble(AC->(double)AC/(double)van).toArray();
				vcb.attribute(this.af.getID(), afs);
				}
			else
				{
				final Double afs[]=new Double[alts.size()];
				Arrays.fill(afs,null);
				vcb.attribute(this.af.getID(),  Arrays.asList(afs));
				}
			}
		}
	
	
    public VCFFamilies()
    	{
    	}
	
    @Override
    protected Logger getLogger() {
    	return LOG;
    	}
    
	@Override
	public int doVcfToVcf(final String inputName, VCFIterator r, final VariantContextWriter w)
		{
		final VCFHeader header= r.getHeader();
		final Pedigree pedigree;
		try {
			pedigree=new PedigreeParser().parse(this.pedigreeFile);
			}
		catch(final Throwable err)
			{
			throw new RuntimeIOException(err);
			}
		
		if(pedigree==null || pedigree.isEmpty())
			{
			throw new RuntimeIOException("Pedigree null/empty");
			}
		final VCFHeader h2=new VCFHeader(header);
		final Map<String,FamilyInfo> famidToFamilyInfo = new HashMap<>();
		
		pedigree.getSamplesInVcfHeader(header).
				forEach(P->{
				final Family pedFamily = P.getFamily();
				FamilyInfo finfo = famidToFamilyInfo.get(pedFamily.getId());
				if(finfo==null) {
					finfo = new FamilyInfo(pedFamily);
					famidToFamilyInfo.put(pedFamily.getId(),finfo);
					}
				finfo.samples.add(P.getId());
			});
		
		famidToFamilyInfo.values().stream().flatMap(F->F.getMetaDataLines().stream()).
			forEach(H->h2.addMetaDataLine(H));
		JVarkitVersion.getInstance().addMetaData(this, h2);
		
		
		w.writeHeader(h2);
		while(r.hasNext())
			{
			final VariantContext ctx = r.next();
			final VariantContextBuilder vcb= new VariantContextBuilder(ctx);
			final List<Allele> alts = ctx.getAlternateAlleles();
			famidToFamilyInfo.values().forEach(F->F.visit(vcb,ctx,alts));
			w.add(vcb.make());	
			}
		w.close();
		return 0;
		}

	
	public static void main(final String[] args)
		{
		new VCFFamilies().instanceMainWithExit(args);
		}

	}
