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
package com.github.lindenb.jvarkit.tools.vcftrios;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BiPredicate;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlType;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.JexlGenotypePredicate;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
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
		generate_doc=false
		)
public class VCFFamilies
	extends Launcher
	{
	private static final  Logger LOG = Logger.build(VCFFamilies.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();

	@XmlType(name="vcffamilies")
	@XmlRootElement(name="vcffamilies")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class CtxWriterFactory 
	implements VariantContextWriterFactory
		{
		private class CtxWriter extends DelegateVariantContextWriter
			{
			private final File pedigreeFile = CtxWriterFactory.this.pedigreeFile;
			private final Map<String,FamilyInfo> famidToFamilyInfo=new HashMap<>();
			private final BiPredicate<VariantContext,Genotype> genotypeFilter = CtxWriterFactory.this.genotypeFilter;

			private class FamilyInfo
				{
				final Pedigree.Family pedFamily;
				private final Set<String> samples=new HashSet<>();
				private final VCFInfoHeaderLine ac;
				private final VCFInfoHeaderLine an;
				private final VCFInfoHeaderLine af;
				
				
				
				FamilyInfo(final Pedigree.Family pedFamily)
					{
					this.pedFamily=pedFamily;
					final String prefix=
							(
							StringUtil.isBlank(CtxWriterFactory.this.prefix)?
							CtxWriterFactory.this.prefix+"_":""
							)+
						this.pedFamily.getId()
						;
					this.ac = new VCFInfoHeaderLine(
							prefix+"_AC",
							VCFHeaderLineCount.A,
							VCFHeaderLineType.Integer,
							"Number of alleles carrying the variant allele for the family "+this.pedFamily.getId()
							);
					this.an = new VCFInfoHeaderLine(
							prefix+"_AN",
							1,
							VCFHeaderLineType.Integer,
							"Number of alleles for the family "+this.pedFamily.getId()
							);
					this.af = new VCFInfoHeaderLine(
							prefix+"_AF",
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
							filter(G->CtxWriter.this.genotypeFilter.test(ctx, G)).
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
								filter(G->CtxWriter.this.genotypeFilter.test(ctx, G)).
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
			
			
			CtxWriter(final VariantContextWriter delegate) {
				super(delegate);
				}
			 
			@Override
			public void writeHeader(final VCFHeader header) {
				final Pedigree pedigree;
				try {
					if(this.pedigreeFile==null)
						{
						pedigree = Pedigree.newParser().parse(header);
						}
					else
						{
						LOG.info("reading pedigree "+this.pedigreeFile);
						pedigree=Pedigree.newParser().parse(this.pedigreeFile);
						}
					}
				catch(final Exception err)
					{
					LOG.error(err);
					throw new RuntimeIOException(err);
					}
				
				if(pedigree==null || pedigree.isEmpty())
					{
					throw new RuntimeIOException("Pedigree null/empty");
					}
				final VCFHeader h2=new VCFHeader(header);
				this.famidToFamilyInfo.clear();
				
				final Set<String> sampleSet = new HashSet<>(header.getSampleNamesInOrder());
				pedigree.getPersons().stream().
						filter(P->sampleSet.contains(P.getId())).
						forEach(P->{
						final Pedigree.Family pedFamily = P.getFamily();
						FamilyInfo finfo = this.famidToFamilyInfo.get(pedFamily.getId());
						if(finfo==null) {
							finfo = new FamilyInfo(pedFamily);
							this.famidToFamilyInfo.put(pedFamily.getId(),finfo);
							}
						finfo.samples.add(P.getId());
					});
				
				this.famidToFamilyInfo.values().stream().flatMap(F->F.getMetaDataLines().stream()).
					forEach(H->h2.addMetaDataLine(H));
				
				super.writeHeader(h2);				
				}
			
			@Override
			public void add(final VariantContext ctx) {
				final VariantContextBuilder vcb= new VariantContextBuilder(ctx);
				final List<Allele> alts = ctx.getAlternateAlleles();
				this.famidToFamilyInfo.values().forEach(F->F.visit(vcb,ctx,alts));
				super.add(vcb.make());				
				}
			
			@Override
			public void close() {
				super.close();
				}
			}
		
		@XmlElement(name="pedigree")
		@Parameter(names={"-p","--ped","--pedigree"},description="Pedigree file. "+Pedigree.OPT_DESCRIPTION)
		private File pedigreeFile = null;
		@Parameter(names={"-pfx","--prefix"},description="VCF header Attribute prefix")
		private String prefix = "";
		@Parameter(names={"-gf","--genotype-filter"},description=JexlGenotypePredicate.PARAMETER_DESCRIPTION,converter=JexlGenotypePredicate.Converter.class)
		private BiPredicate<VariantContext,Genotype> genotypeFilter = JexlGenotypePredicate.create("!g.isFiltered()");
		
		
		@Override
		public VariantContextWriter open(final VariantContextWriter delegate) {
			return new CtxWriter(delegate);
			}
		}
	
    public VCFFamilies()
    	{
    	}
		
	@Override
	public int doVcfToVcf(final String inputName, VCFIterator r, final VariantContextWriter delegate)
		{
		final VariantContextWriter out = this.component.open(delegate);
		final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(r.getHeader()).logger(LOG);
		out.writeHeader(r.getHeader());
		while(r.hasNext())
			{
			out.add(progress.watch(r.next()));
			}
		out.close();
		progress.finish();
		return 0;
		}

	@Override
	public int doWork(final List<String> args) {
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
	
	public static void main(final String[] args)
		{
		new VCFFamilies().instanceMainWithExit(args);
		}

	}
