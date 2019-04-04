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
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.PostponedVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VariantAttributesRecalculator;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfTools;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
/*
BEGIN_DOC

## SNPEFF/VEP Annotations

this tool will try to split the VEP or the VCF annotation for each allele.


## Example

Exac contains multi-ALT  variants:

```bash
$ gunzip -c ExAC.r0.3.sites.vep.vcf.gz | grep rs3828049

1	889238	rs3828049	G	A,C	8422863.10	PASS	AC=6926,3;AC_AFR=220,0;AC_AMR=485,1;AC_Adj=6890,3;AC_EAS=746,0;AC_FIN=259,0;AC_Het=6442,3,0;AC_Hom=224,0;AC_NFE=3856,0;AC_OTH=41,0;AC_SAS=1283,2;AF=0.057,2.472e-05;AN=121358;AN_AFR=10148;AN_AMR=11522;AN_Adj=119272;AN_EAS=8582;AN_FIN=6358;AN_NFE=65282;AN_OTH=876;AN_SAS=16504;(...)

```

processed with this tools:
```
$ java -jar dist/vcfmulti2oneallele.jar  ExAC.r0.3.sites.vep.vcf.gz   | grep rs3828049

1	889238	rs3828049	G	A	8422863.10	PASS	AC=6926;AC_AFR=220;AC_AMR=485;AC_Adj=6890;AC_EAS=746;AC_FIN=259;AC_Het=6442;AC_Hom=224;AC_NFE=3856;AC_OTH=41;AC_SAS=1283;AF=0.057;AN=121358;AN_AFR=10148;AN_AMR=11522;AN_Adj=119272;AN_EAS=8582;AN_FIN=6358;AN_NFE=65282;AN_OTH=876;AN_SAS=16504;BaseQRankSum=-2.170e-01;VCF_MULTIALLELIC_SRC=A|C;(...)
1	889238	rs3828049	G	C	8422863.10	PASS	AC=3;AC_AFR=0;AC_AMR=1;AC_Adj=3;AC_EAS=0;AC_FIN=0;AC_Het=3;AC_Hom=0;AC_NFE=0;AC_OTH=0;AC_SAS=2;AF=2.472e-05;AN=121358;AN_AFR=10148;AN_AMR=11522;AN_Adj=119272;AN_EAS=8582;AN_FIN=6358;AN_NFE=65282;AN_OTH=876;AN_SAS=16504;VCF_MULTIALLELIC_SRC=A|C;(....)
```

## History

* 20170606 added support for VCFHeaderLineCount.R

END_DOC
 */
@Program(
		name="vcfmulti2oneallele",
		description="'one variant with N ALT alleles' to 'N variants with one ALT'",
		keywords={"vcf"}
		)
public class VcfMultiToOneAllele
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfMultiToOneAllele.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@ParametersDelegate
	private PostponedVariantContextWriter.WritingVcfConfig writingVcfArgs = new PostponedVariantContextWriter.WritingVcfConfig();
	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();

	
	public static class CtxWriterFactory 
		implements VariantContextWriterFactory
		{
		@Parameter(names={"-p","--samples"},description="print sample genotypes.")
		private boolean print_samples = false;
		@Parameter(names={"-r","--rmAtt"},description="[20161110]: after merging with GATK CombineVariants there can have problemes with INFO/type='A' present in vcf1 but not in vcf2, and multiallelelic variants. This option delete the attributes having such problems.")
		private boolean rmErrorAttributes = false;
		@Parameter(names={"-highest","--highest"},description="[20170723]: Use  Allele With Highest Allele Count, discard/replace the other")
		private boolean useAltAlleleWithHighestAlleleCount = false;
		@Parameter(names={"-tag","--tag"},description="Info field name that will be added to recall the original alleles.")
		private String TAG="VCF_MULTIALLELIC_SRC";
		@Parameter(names={"--ignoreMissingInfoDecl"},description="Ignore error when a variant INFO is missing a definition in the VCF header.")
		private boolean ignoreMissingInfoDecl=false;
		@Parameter(names={"--addNoVariant"},description="Print Variants without ALT allele")
		private boolean addNoVariant=false;
		@Parameter(names={"--skipSpanningDeletions"},description="Skip Alt Spanning deletion alleles "+Allele.SPAN_DEL_STRING)
		private boolean skipSpanningDeletion=false;
		@Parameter(names={"--replaceWith"},description="When replacing an alternative allele, replace it with REF or current ALT allele.")
		private ReplaceWith replaceWith=ReplaceWith.REF;
		@Parameter(names={"--disableHomVarAlt"},description="by default is a genotype is homvar for an external ALT ('2/2'), it will be set to ./. (no call). Setting this option will replace the current allele.")
		private boolean disableHomVarAlt=true;
		@ParametersDelegate
		private VariantAttributesRecalculator recalculator = new VariantAttributesRecalculator();

		
		private class CtxWriter extends DelegateVariantContextWriter
			{
			private VcfTools tools;
			private VCFHeader header2 = null;
			private final Function<List<Allele>, String> altListToString = alternateAlleles -> alternateAlleles.stream().
					map(A->A.getDisplayString()).
					collect(Collectors.joining("|"))
					;
			private List<String> sample_names = Collections.emptyList();
			
			CtxWriter(final VariantContextWriter delegate) {
				super(delegate);
				}
			
			@Override
			public void writeHeader(VCFHeader header) {
				this.tools = new VcfTools(header);
				this.sample_names=header.getSampleNamesInOrder();
				final Set<VCFHeaderLine> metaData=new HashSet<>(header.getMetaDataInInputOrder());
				//addMetaData(metaData);		
				metaData.add(new VCFInfoHeaderLine(
						CtxWriterFactory.this.TAG,
						1,
						VCFHeaderLineType.String,
						"The variant was processed with VcfMultiAlleleToOneAllele and contained the following alleles."));
				
				if(!CtxWriterFactory.this.print_samples)
					{
					this.sample_names = Collections.emptyList();
					this.header2 = new VCFHeader(
							metaData,
							this.sample_names
							);
					}
				else
					{
					this.header2 = new VCFHeader(
							metaData,
							sample_names
							);
					}	
				CtxWriterFactory.this.recalculator.setHeader(this.header2);
				super.writeHeader(this.header2);
				}
			
			@Override
			public void add(final VariantContext ctx) {
				final List<Allele> alternateAlleles = new ArrayList<>(ctx.getAlternateAlleles());
				if(alternateAlleles.isEmpty())
					{
					if(CtxWriterFactory.this.addNoVariant)
						{
						if(!CtxWriterFactory.this.print_samples)
							{
							final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
							vcb.noGenotypes();
							super.add(vcb.make());
							}
						else
							{
							super.add(ctx);
							}						
						}
					else
						{
						LOG.warn("Remove no ALT variant:"+ctx);
						}
					}
				else if(alternateAlleles.size()==1)
					{		
					if(CtxWriterFactory.this.skipSpanningDeletion && alternateAlleles.get(0).equals(Allele.SPAN_DEL))
						{
						return;
						}
					
					if(!CtxWriterFactory.this.print_samples)
						{
						final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
						vcb.noGenotypes();
						super.add(vcb.make());
						}
					else
						{
						super.add(ctx);
						}
					}
				else
					{
					final Allele highest = (CtxWriterFactory.this.useAltAlleleWithHighestAlleleCount?
							ctx.getAltAlleleWithHighestAlleleCount():
							null
							);

					
					for(int alternateIndex=0;alternateIndex< alternateAlleles.size();++alternateIndex)
						{
						final Allele the_allele = alternateAlleles.get(alternateIndex);
						if(highest!=null && !highest.equals(the_allele)) continue;
						if(CtxWriterFactory.this.skipSpanningDeletion && the_allele.equals(Allele.SPAN_DEL))
							{
							continue;
							}
						
						
						final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
						vcb.alleles(Arrays.asList(ctx.getReference(),the_allele));
						vcb.attributes(makeAttributes(this.header2,this.tools,ctx,the_allele));						
						vcb.attribute(CtxWriterFactory.this.TAG,altListToString.apply(alternateAlleles));
						
						if(!print_samples)
							{
							vcb.noGenotypes();
							}
						else
							{							
							final Allele replaceAlleleWith;
							switch(CtxWriterFactory.this.replaceWith)
								{
								case ALT: replaceAlleleWith= the_allele;break;
								case REF: replaceAlleleWith= ctx.getReference();break;
								default: throw new IllegalStateException();
								}
							vcb.genotypes(makeGenotypes(ctx, sample_names, the_allele, replaceAlleleWith));
							}
						super.add(CtxWriterFactory.this.recalculator.apply(vcb.make()));
						}
					}
				
				}
			
			
			
			}
		
		@Override
		public VariantContextWriter open(VariantContextWriter delegate) {
			return new CtxWriter(delegate);
			}
		
		
		 private Map<String,Object> makeAttributes(
				 final VCFHeader header,
				 final VcfTools vcfTools,
				 final VariantContext ctx,
				 final Allele the_allele
				 )
		 	{
			final int alleleIndex = ctx.getAlleleIndex(the_allele);
			 
			if(alleleIndex==0) throw new IllegalStateException("alleleIndex==0 (REF)");
			if(alleleIndex==-1) throw new IllegalStateException("alleleIndex<0");
			final Map<String,Object> attributes = ctx.getAttributes();
			final Map<String,Object> newAttributes = new HashMap<>(attributes.size());
				

				
			for(final String attid:attributes.keySet())
				{
				final VCFInfoHeaderLine info = header.getInfoHeaderLine(attid);
				if(info==null) {
					final String msg = "Cannot get header INFO tag="+attid+" at "+ctx.getContig()+":"+ctx.getStart();
					if(!this.ignoreMissingInfoDecl)
						{
						throw new JvarkitException.FileFormatError(msg);
						}
					else
						{
						LOG.warning(msg);
						continue;
						}
					}
				// get ANN specific annotations
				if(info.getID().equals(vcfTools.getAnnPredictionParser().getTag()))
					{
					final List<String> L=
						vcfTools.getAnnPredictionParser().getPredictions(ctx).
						stream().
						filter(P->the_allele.getDisplayString().equals(P.getAllele())).
						map(P->P.getOriginalAttributeAsString()).
						collect(Collectors.toList());
					if(L.isEmpty())
						{
						if(!the_allele.equals(Allele.SPAN_DEL)) {
							LOG.warning("for "+ctx.getContig()+":"+ctx.getStart()+". No ANN Prediction for "+the_allele +" in "+
									ctx.getAttributeAsList(vcfTools.getVepPredictionParser().getTag()));
							}
						}
					else
						{
						newAttributes.put(info.getID(), L);
						}
					continue;
					}
				// get VEP specific annotations
				if(info.getID().equals(vcfTools.getVepPredictionParser().getTag()))
					{
					List<String> L=
							vcfTools.getVepPredictionParser().getPredictions(ctx).
							stream().
							filter(P->the_allele.getDisplayString().equals(P.getAlleleStr())).
							map(P->P.getOriginalAttributeAsString()).
							collect(Collectors.toList());
					
					if(L.isEmpty() && ctx.isIndel())
						{
						L=  vcfTools.getVepPredictionParser().getPredictions(ctx).
							stream().
							filter(P->P.getAlleleStr().equals(VepPredictionParser.INDEL_SYMBOL_STR)).
							map(P->P.getOriginalAttributeAsString()).
							collect(Collectors.toList());
						}
					
					if(L.isEmpty())
						{
						if(!the_allele.equals(Allele.SPAN_DEL)) {
							LOG.warning("for "+ctx.getContig()+":"+ctx.getStart()+". No Vep Prediction for "+the_allele +" in "+
									ctx.getAttributeAsList(vcfTools.getVepPredictionParser().getTag()));
							}
						}
					else
						{
						newAttributes.put(info.getID(), L);
						}
					continue;
					}
				final VCFHeaderLineCount lineCount = info.getCountType();

				if(lineCount!=VCFHeaderLineCount.A && lineCount!=VCFHeaderLineCount.R)
					{
					newAttributes.put(attid, attributes.get(attid));
					continue;
					}
				final Object o = 	attributes.get(attid);
				
				if(!(o instanceof List)) {
					final String msg="For INFO tag="+attid+" got "+o.getClass()+" instead of List in "+ctx;
					if(this.rmErrorAttributes)
						{
						LOG.warn("remove this attribute : "+msg);
						continue;
						}
					else
						{
						throw new JvarkitException.FileFormatError(msg);
						}				
					}
				@SuppressWarnings("rawtypes")
				final List list = (List)o;
				if(ctx.getNAlleles() != list.size()+  ( lineCount.equals(VCFHeaderLineCount.A) ? 1 : 0 ) ) 
					{
					final String msg= ctx.getContig()+":"+ctx.getStart()+" : For INFO tag="+attid+" got "+ctx.getNAlleles()+" ALLELES, incompatible with "+list.toString();
					if(this.rmErrorAttributes)
						{
						LOG.warn("remove this attribute : "+msg);
						continue;
						}
					else
						{
						throw new JvarkitException.FileFormatError(msg);
						}
					}
				else if(lineCount.equals(VCFHeaderLineCount.R))
					{
					newAttributes.put(attid, Arrays.asList(
							list.get(0)/* REF */,
							list.get(alleleIndex))
							);	
					}
				else // VCFHeaderLineCount.A
					{	
					newAttributes.put(attid, list.get(alleleIndex-1) /* -1 because the index is in the total allele list, including REF=0 */);	
					}
				}
			return	newAttributes;
		 	}
		 
		 
		 private List<Genotype> makeGenotypes(
				 final VariantContext ctx,
				 final List<String> sample_names,
				 final Allele theAllele,
				 final Allele replaceWith
				 )
		 		{
				final List<Genotype> genotypes=new ArrayList<>(sample_names.size());
				
				for(final String sampleName: sample_names)
					{							
					final Genotype g= ctx.getGenotype(sampleName);
					
					if( !disableHomVarAlt &&
						g.isCalled() && 
						!g.getAlleles().stream().
						filter(A->!(A.isNoCall() || A.isReference() || A.equals(theAllele) || A.equals(replaceWith))).
						collect(Collectors.toSet()).
						isEmpty() // only contains the 'other alleles'
						)
						{
						genotypes.add(GenotypeBuilder.createMissing(sampleName, g.getPloidy()));
						continue;
						}
					
					final List<Allele> newAllelesList = g.getAlleles().stream().
							map(A->(A.isNoCall() || A.isReference() || A.equals(theAllele)?A:replaceWith)).
							collect(Collectors.toList())
							;
					final GenotypeBuilder gb =new GenotypeBuilder(
							g.getSampleName(),
							newAllelesList
							);
					if(g.hasDP()) gb.DP(g.getDP());
					if(g.hasGQ()) gb.GQ(g.getGQ());
					if(g.isFiltered()) gb.filter(g.getFilters());
					if(g.hasAD())
						{
						final int srcads[]=g.getAD();
						if(srcads==null || srcads.length==0 ||  newAllelesList.size()!=2) 
							{
							// nothing
							}
						else if(newAllelesList.get(0).equals(newAllelesList.get(1)) )
							{
							//HOM_REF
							if(newAllelesList.get(0).isReference())
								{
								gb.AD(new int[] {srcads[0],0});
								gb.DP(srcads[0]);
								}
							//HOM_VAR
							else if(newAllelesList.get(0).equals(replaceWith))
								{
								final int aidx = ctx.getAlleleIndex(replaceWith);
								if(aidx>=0 && aidx < srcads.length)
									{
									gb.AD(new int[] {0,srcads[aidx]});
									gb.DP(srcads[aidx]);
									}
								}
							}
						//HET
						else if(
							(newAllelesList.get(0).isReference() && newAllelesList.get(1).equals(replaceWith)) ||
							(newAllelesList.get(1).isReference() && newAllelesList.get(0).equals(replaceWith)))
							{
							final int aidx = ctx.getAlleleIndex(replaceWith);
							if(aidx>=0 && aidx < srcads.length)
								{
								gb.AD(new int[] {srcads[0],srcads[aidx]});
								gb.DP(srcads[0]+srcads[aidx]);
								}
							}
							
						}

					genotypes.add(gb.make());
					}
			return genotypes;
		 	}
		}
	
	enum ReplaceWith {REF,ALT};

	 public VcfMultiToOneAllele()
		{
		}

	@Override
	public int doVcfToVcf(
			final String inputName,
			final VCFIterator r,
			final VariantContextWriter delegate
			)   
		{
		final VariantContextWriter w = this.component.open(delegate);
		w.writeHeader(r.getHeader());
		final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(r.getHeader()).logger(LOG);
		while(r.hasNext())
			{
			w.add(progress.watch(r.next()));
			}
		progress.finish();
		w.close();
		return 0;
		}

	@Override
	public int doWork(final List<String> args) {
		try 
			{
			if(this.component.initialize()!=0) {
				return -1;
				}
			return doVcfToVcf(args,this.outputFile);
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.component);
			}
		}
	
	@Override
	protected VariantContextWriter openVariantContextWriter(final File outorNull) throws IOException {
		return new PostponedVariantContextWriter(this.writingVcfArgs,stdout(),this.outputFile);
		}

	
	public static void main(final String[] args)
		{
		new VcfMultiToOneAllele().instanceMainWithExit(args);
		}
	
	}
