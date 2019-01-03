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
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlType;

import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.PostponedVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfTools;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;


/**
BEGIN_DOC

## Motivation

when using gatk SelectVariants with sample names (-sn) some alleles specific of the samples than have been removed, remain in the vcf.

## SNPEFF / VEP

this tool removes unused annotations from SNPEFF(ANN=) and VEP.

## Example

```bash
$ cat in.vcf
(...)
chr1	7358	.	ACTT	*,A	1313.61	PASS	AC=0,10;AF=0,0.005828;AN=1716

$ java -jar dist/vcfremoveunusedalt.jar  in.vcf | grep -w 17358 -m1
chr1	7358	.	ACTT	A	1313.61	PASS	AC=10;AF=0.005828;AN=1716
```
END_DOC

 */
@SuppressWarnings("unused")
@Program(name="vcfremoveunusedalt",
	description="Remove unused ALT allele if there is no genotype with this alt, or there is no sample but AC=0",
	keywords={"vcf","genotype"}
)
public class VcfRemoveUnusedAlt extends Launcher {
	private static final Logger LOG=Logger.build(VcfRemoveUnusedAlt.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT,required=false)
	private File output=null;
	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();
	
	
	@XmlType(name="vcfremoveunusedalt")
	@XmlRootElement(name="vcfremoveunusedalt")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class CtxWriterFactory 
		implements VariantContextWriterFactory
			{
			@Parameter(names={"-onespan","--onespan"},description="Don't print the variant if the only remaining allele is  '"+Allele.SPAN_DEL_STRING+"'")
			private boolean no_span_allele =false;
			@Parameter(names={"-neverspan","--neverspan"},description="Remove ALL spanning deletions '"+Allele.SPAN_DEL_STRING+"'. VCF must have no genotype.")
			private boolean never_span_allele =false;
		
			private class CtxWriter extends DelegateVariantContextWriter
				{
				private VcfTools vcfTools = null;
				private VCFHeader header;
				private final boolean no_span_allele = CtxWriterFactory.this.no_span_allele;
				private final boolean never_span_allele = CtxWriterFactory.this.never_span_allele;
				private long nChanges = 0L;
				CtxWriter(final VariantContextWriter delegate) {
					super(delegate);
					}
				@Override
				public void writeHeader(final VCFHeader header) {
					this.header=header;
					this.vcfTools = new VcfTools(header);
					
					if(this.never_span_allele && header.getNGenotypeSamples()>0)
						{
						throw new JvarkitException.UserError("cannot remove spanning allele where there are some genotypes");
						}
					
					super.writeHeader(header);
					}
				
				@Override
				public void add(final VariantContext ctx) {
					final List<Allele> alts = ctx.getAlternateAlleles();
					if((this.no_span_allele || this.never_span_allele) && alts.size()==1 && alts.get(0).equals(Allele.SPAN_DEL)) {
						return;
						}
					
					
					
					final Set<Allele> toRemove= new HashSet<>();
					if(ctx.getNSamples()==0) {
						
						if(this.never_span_allele && alts.contains(Allele.SPAN_DEL))
							{
							toRemove.add(Allele.SPAN_DEL);
							}
						
						if(ctx.hasAttribute(VCFConstants.ALLELE_COUNT_KEY))
							{
							final List<Object> acStrList = ctx.getAttributeAsList(VCFConstants.ALLELE_COUNT_KEY);
							if(acStrList.size()!=alts.size()) {
								LOG.error("Illegal number of key "+VCFConstants.ALLELE_COUNT_KEY+" at "+ctx);
								return;
								}
							for(int alt_idx=0;alt_idx< alts.size();++alt_idx)
								{
								final Object o=acStrList.get(alt_idx);
								final String ostr = String.valueOf(o);
								if(o==null || ostr.equals(".") || ostr.equals("0") || Integer.parseInt(ostr)==0)
									{
									toRemove.add(alts.get(alt_idx));
									}
								}
							}
						}
					else
						{
						for(final Allele alt:alts) {
							if(ctx.getGenotypes().stream().
									anyMatch(G->G.getAlleles().contains(alt))) continue;
							toRemove.add(alt);
							}
						}
					
					if(toRemove.size()==alts.size()) {
						return;
						}
					
					if(toRemove.isEmpty()) {
						super.add(ctx);
						return;
						}
					
					final Set<Integer> indices= toRemove.stream().
							map(A->ctx.getAlleleIndex(A)).
							collect(Collectors.toSet());
					if(indices.contains(0)) throw new IllegalStateException();
					if(indices.contains(-1)) throw new IllegalStateException();
					
					
					final VariantContextBuilder gb=new VariantContextBuilder(ctx);
					final List<Allele> newAlleles = ctx.getAlleles().
							stream().
							filter(A->!toRemove.contains(A)).
							collect(Collectors.toList())
							;
					if(this.no_span_allele && alts.size()==2 && alts.get(1).equals(Allele.SPAN_DEL)) {
						return;
						}
					gb.alleles(newAlleles);
					
					// SnpEff AN
					if(ctx.hasAttribute(this.vcfTools.getAnnPredictionParser().getTag()))
						{
						final Set<String> toRemoveStr = toRemove.stream().map(S->S.getDisplayString()).collect(Collectors.toSet());
						final List<String> newpreds = 
								this.vcfTools.getAnnPredictionParser().
								getPredictions(ctx).stream().
									filter(P->!toRemoveStr.contains(P.getAllele())).
									map(P->P.getOriginalAttributeAsString()).
									collect(Collectors.toList());
						
						gb.attribute(this.vcfTools.getAnnPredictionParser().getTag(),newpreds);
						}
					
					// VEP
					if(ctx.hasAttribute(this.vcfTools.getVepPredictionParser().getTag()))
						{
						final Set<String> toRemoveStr = toRemove.stream().map(S->S.getDisplayString()).collect(Collectors.toSet());
						final List<String> newpreds = 
								this.vcfTools.getVepPredictionParser().
								getPredictions(ctx).stream().
									filter(P->!toRemoveStr.contains(P.getAlleleStr())).
									map(P->P.getOriginalAttributeAsString()).
									collect(Collectors.toList());
						
						gb.attribute(this.vcfTools.getVepPredictionParser().getTag(),newpreds);
						}
					
					for(final VCFInfoHeaderLine vcfInfoHeaderLine:this.header.getInfoHeaderLines()) {
						switch(vcfInfoHeaderLine.getCountType())
							{	
							case A: case R: break;
							default: continue;
							}
						
						if(!ctx.hasAttribute(vcfInfoHeaderLine.getID())) continue;
						final List<Object> att = new ArrayList<>(ctx.getAttributeAsList(vcfInfoHeaderLine.getID()));
						final List<Object> newatt = new ArrayList<>(att.size()); 

						if(vcfInfoHeaderLine.getCountType()==VCFHeaderLineCount.R)
							{
							newatt.add(att.get(0));
							for(int x= 1;x < att.size();++x)
								{
								if(indices.contains(x)) continue;
								newatt.add(att.get(x));
								}
							}
						else if(vcfInfoHeaderLine.getCountType()==VCFHeaderLineCount.A)
							{
							for(int x= 0;x < att.size();++x)
								{
								if(indices.contains(x+1)) continue;
								newatt.add(att.get(x));
								}
							}
						else
							{
							throw new IllegalStateException();
							}
						gb.attribute(vcfInfoHeaderLine.getID(),newatt);
						}
					this.nChanges++;
					super.add(gb.make());
					}
				@Override
				public void close() {
					LOG.info("number of changes "+this.nChanges);
					super.close();
					}
				}

			@Override
			public VariantContextWriter open(final VariantContextWriter delegate) {
				return new CtxWriter(delegate);
				}
			}
		
		public VcfRemoveUnusedAlt()
			{
			}
		
		@Override 
		protected int doVcfToVcf(
				final String inputName, 
				final VcfIterator in,
				final VariantContextWriter delegate) {
			try
				{
				final VariantContextWriter out  = this.component.open(delegate);
				out.writeHeader(in.getHeader());
				final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(in.getHeader()).logger(LOG);
				while(in.hasNext())
					{
					out.add(progress.watch(in.next()));
					}
				progress.finish();
				out.close();
				return 0;
				}
			catch(final Exception err)
				{
				LOG.error(err);
				return -1;
				}
			}
		
		@Override
		public int doWork(final List<String> args) {
			return doVcfToVcf(args,this.output);
			}
	
		public static void main(final String[] args)
			{
			new VcfRemoveUnusedAlt().instanceMainWithExit(args);
			}
		}
