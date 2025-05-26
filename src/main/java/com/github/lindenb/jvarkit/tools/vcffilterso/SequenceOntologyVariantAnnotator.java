package com.github.lindenb.jvarkit.tools.vcffilterso;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.vcf.VariantAttributesRecalculator;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.BcfToolsPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.BcfToolsPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;
import com.github.lindenb.jvarkit.variant.VariantAnnotator;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFIterator;

public class SequenceOntologyVariantAnnotator implements VariantAnnotator {
	private static final Logger LOG = Logger.of(VcfFilterSequenceOntology.class);
	private static final String GT_FILTER_RESET_TO_NOCALL="NO_CALL";
	
	private boolean invertSoTerms = false;

	private boolean disableReasoning = false;

	private String filterIn = null;

	private String filterOut = null;

	private boolean removeUnusedAttribute = false;

	private boolean removeVariantIfNoMoreAttribute = false;

	private String filterGenotypesStr=null;
	
	private SequenceOntologyTree _sequenceOntologyTree = null;

	/* all sequence terms */
	private final Set<SequenceOntologyTree.Term> user_terms=new HashSet<>();

	
	private final List<AbstractPredictionHandler> predictionHandlers = new ArrayList<>();				

	
	public void setSequenceOntologyTree(SequenceOntologyTree sequenceOntologyTree) {
		this._sequenceOntologyTree = sequenceOntologyTree;
		}
	
	public SequenceOntologyTree getSequenceOntologyTree() {
		if(_sequenceOntologyTree==null) {
			_sequenceOntologyTree = SequenceOntologyTree.createDefault();
			}
		return _sequenceOntologyTree;
		}
	
	public boolean isRemoveUnusedAttribute() {
		return removeUnusedAttribute;
		}

	private abstract class AbstractPredictionHandler
		{
		final List<Object> predStrings = new ArrayList<>();
		final Set<Allele> matching_alleles = new HashSet<>();
		boolean keepFlag = false;
		abstract boolean isValid();
		abstract String getTag();
		
		AbstractPredictionHandler reset(final VariantContext ctx) {
			this.predStrings.clear();
			this.keepFlag = false;
			this.matching_alleles.clear();
			this.matching_alleles.add(ctx.getReference());
			return this;
			}
		boolean supportFilterAlleles() { return false;}
		abstract AbstractPredictionHandler visit(final VariantContext ctx,final VariantContextBuilder vcb);
		void updateInfo(final VariantContext ctx,final VariantContextBuilder vcb)
			{
			if(isRemoveUnusedAttribute()) {
				vcb.rmAttribute(this.getTag());
				if(!this.predStrings.isEmpty()) {
					vcb.attribute(this.getTag(), this.predStrings);
					}
				}
			}	
		}
	
	/** handler for VEP */
	private class VepPredictionHandler extends AbstractPredictionHandler
		{
		private final VepPredictionParser parser;
		VepPredictionHandler(final VCFHeader header)
			{
			this.parser =new VepPredictionParserFactory().
					header(header).get().
					sequenceOntologyTree(getSequenceOntologyTree());
			}
		@Override
		String getTag() {
			return this.parser.getTag();
			}
		@Override
		boolean isValid() { return this.parser.isValid();}
		@Override
		boolean supportFilterAlleles() {
			return true;
			}
		AbstractPredictionHandler visit(final VariantContext ctx,final VariantContextBuilder vcb)
			{
			if(!ctx.hasAttribute(this.getTag())) return this;
			for(final VepPredictionParser.VepPrediction pred : this.parser.getPredictions(ctx))
				{

				if(pred==null) continue;
				if(hasUserTemLabel(pred.getSOTerms()))
					{
					if(isRecodingGenotypes()) {
						if(pred.getAllele()!=null) this.matching_alleles.add(pred.getAllele());
						}
					this.keepFlag=true;
					if(isRemoveUnusedAttribute()) {
						this.predStrings.add(pred.getOriginalAttributeAsString());
						}
					}
				else
					{
					//nothing
					}
				}
			updateInfo(ctx,vcb);
			return this;
			}
		}
				
	/** handler for SNPEFF */
	private class SnpEffPredictionHandler extends AbstractPredictionHandler
		{
		private final SnpEffPredictionParser parser;
		SnpEffPredictionHandler(final VCFHeader header)
			{
			this.parser =new SnpEffPredictionParserFactory().
					header(header).get().
					sequenceOntologyTree(getSequenceOntologyTree());
			}
		
		@Override
		boolean supportFilterAlleles() {
			return true;
			}
		@Override
		String getTag() {
			return this.parser.getTag();
			}
		@Override
		boolean isValid() { return this.parser.isValid();}
		@Override
		AbstractPredictionHandler visit(final VariantContext ctx,final VariantContextBuilder vcb)
			{
			if(!ctx.hasAttribute(this.getTag())) return this;
			for(final SnpEffPredictionParser.SnpEffPrediction pred : this.parser.getPredictions(ctx))
				{
				if(pred==null) continue;
				if(hasUserTemLabel(pred.getSOTerms()))
					{
					if(isRecodingGenotypes()) {
						final Allele alt = pred.getAllele();
						if(alt!=null) this.matching_alleles.add(alt);
						}
					this.keepFlag=true;
					if(isRemoveUnusedAttribute()) {
						this.predStrings.add(pred.getOriginalAttributeAsString());
						}
					}
				}
			updateInfo(ctx,vcb);
			return this;
			}
		}
				
	/** prediction for ANN  */
	private class AnnPredictionHandler extends AbstractPredictionHandler
		{
		private final AnnPredictionParser parser;
		AnnPredictionHandler(final VCFHeader header)
			{
			this.parser = new AnnPredictionParserFactory().
					header(header).
					get().
					sequenceOntologyTree(getSequenceOntologyTree());
			}
		@Override
		String getTag() {
			return this.parser.getTag();
			}
		@Override
		boolean isValid() { return this.parser.isValid();}
		@Override
		boolean supportFilterAlleles() {
			return true;
			}
		
		AbstractPredictionHandler visit(final VariantContext ctx,final VariantContextBuilder vcb)
			{
			if(!ctx.hasAttribute(this.getTag())) return this;
			for(final AnnPredictionParser.AnnPrediction pred : this.parser.getPredictions(ctx))
				{
				if(pred==null) continue;
				if(hasUserTemLabel(pred.getSOTerms()))
					{
					if(isRecodingGenotypes() && !StringUtil.isBlank(pred.getAllele()))
						{
						this.matching_alleles.add(Allele.create(pred.getAllele(),false));
						}
					this.keepFlag=true;
					if(isRemoveUnusedAttribute()) {
						this.predStrings.add(pred.getOriginalAttributeAsString());
						}
					}
				}
			updateInfo(ctx,vcb);
			return this;
			}
		}
				
	
	/** prediction for BCFTools CSQ */
	private class BcftoolsCsqPredictionHandler extends AbstractPredictionHandler
		{
		private final BcfToolsPredictionParser parser;
		BcftoolsCsqPredictionHandler(final VCFHeader header)
			{
			this.parser = new BcfToolsPredictionParserFactory().header(header).get().
					sequenceOntologyTree(getSequenceOntologyTree());
			}
		@Override
		String getTag() {
			return this.parser.getTag();
			}
		@Override
		boolean isValid() { return this.parser.isValid();}
		@Override
		boolean supportFilterAlleles() {
			return true;
			}
		
		AbstractPredictionHandler visit(final VariantContext ctx,final VariantContextBuilder vcb)
			{
			if(!ctx.hasAttribute(this.getTag())) return this;
			for(final BcfToolsPredictionParser.BcfToolsPrediction pred : this.parser.getPredictions(ctx))
				{
				if(pred==null) continue;
				if(hasUserTemLabel(pred.getSOTerms()))
					{
					if(isRecodingGenotypes() && !StringUtil.isBlank(pred.getAllele()))
						{
						this.matching_alleles.add(Allele.create(pred.getAllele(),false));
						}
					this.keepFlag=true;
					if(isRemoveUnusedAttribute()) {
						this.predStrings.add(pred.getOriginalAttributeAsString());
						}
					}
				}
			updateInfo(ctx,vcb);
			return this;
			}
		}
	
	
		protected int doVcfToVcf(String inputName, VCFIterator iter, VariantContextWriter out) {
			final List<AbstractPredictionHandler> predictionHandlers = new ArrayList<>();				
			
			final VCFHeader header0 = iter.getHeader();
			
			
			final VCFHeader header2= new VCFHeader(header0);
			
			final String termlist = String.join(", ",this.user_terms.stream().
						map(S->S.getAcn()+"("+S.getLabel()+")").
						collect(Collectors.toSet()))
						;
			if(!StringUtils.isBlank(this.filterIn)) {
				header2.addMetaDataLine(new VCFFilterHeaderLine(this.filterIn,
						"Variant having SO terms:"+ termlist));
				}
			if(!StringUtils.isBlank(this.filterOut)) {
				header2.addMetaDataLine(new VCFFilterHeaderLine(this.filterOut,
						"Variant non having SO terms :" + termlist));
				}
			
			if(!StringUtil.isBlank(this.filterGenotypesStr) &&
				!GT_FILTER_RESET_TO_NOCALL.equals(this.filterGenotypesStr))
				{
				header2.addMetaDataLine(new VCFFormatHeaderLine(
						VCFConstants.GENOTYPE_FILTER_KEY,
						1,VCFHeaderLineType.String,
						"Genotype was filterered by vcffilterso : "+this.filterGenotypesStr
						));
				}
			
			AbstractPredictionHandler ph = new VepPredictionHandler(header0);
			if(ph.isValid()) predictionHandlers.add(ph);
			ph = new SnpEffPredictionHandler(header0);
			if(ph.isValid()) predictionHandlers.add(ph);
			ph = new AnnPredictionHandler(header0);
			if(ph.isValid()) predictionHandlers.add(ph);
			ph = new BcftoolsCsqPredictionHandler(header0);
			if(ph.isValid()) predictionHandlers.add(ph);
			
			this.recalculator.setHeader(header2);
			out.writeHeader(header2);
			}
			
	private VariantContext annotate1(final VariantContext ctx) {
		final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
		
		/* loop over each handler to detect the matching predictions */
		for(final AbstractPredictionHandler handler: predictionHandlers) {
			handler.reset(ctx).visit(ctx,vcb);
			}
		
		/* FILTER genotypes having NO causal ATL allele */
		if(isRecodingGenotypes() && !predictionHandlers.isEmpty()) {
			/* samples to be FILTERED */
			final Set<String> invalidSamples = new HashSet<>(ctx.getNSamples());
			/* loop over each genotype */
			for(final String sample: ctx.getSampleNames())
				{
				boolean sample_is_ok = false;
				for(final AbstractPredictionHandler handler: predictionHandlers) {							
					if(!handler.supportFilterAlleles()) continue;
					final Genotype gt = ctx.getGenotype(sample);
					if(gt==null || gt.isNoCall() || gt.isHomRef() || gt.isFiltered())
						{
						sample_is_ok = true;
						break;
						}
					if(gt.getAlleles().stream().
						filter(A->A.isCalled() && !A.isReference()).
						anyMatch(A->handler.matching_alleles.contains(A)))
						{
						sample_is_ok = true;
						break;
						}
					}
				if(!sample_is_ok) {
					invalidSamples.add(sample);
					}
				}
			final Function<Genotype,Genotype> convertGt = G->{
				/* sample is not invalid */
				if(!invalidSamples.contains(G.getSampleName()))
					{
					return G;
					}
				/* sample is invalid and we reset to NO_CALL (./.) */
				else if(GT_FILTER_RESET_TO_NOCALL.equals(this.filterGenotypesStr))
					{
					return GenotypeBuilder.createMissing(G.getSampleName(), G.getPloidy());
					}
				/* sample is invalid and we set the GT FILTER */
				else
					{
					
					return new GenotypeBuilder(G).filter(this.filterGenotypesStr).make();
					}	
				};
			
			/* update genotypes */
			vcb.genotypes(ctx.getGenotypes().stream().
					map(convertGt).
					collect(Collectors.toList()));
			}
		
		boolean variant_has_one_matching_pred = true;
		if(predictionHandlers.isEmpty() && !this.removeUnusedAttribute) {
			variant_has_one_matching_pred = false;
			}
		

		/* all attributes have been removed ? should we keep this variant ?*/
		if( this.removeUnusedAttribute &&
			predictionHandlers.stream().allMatch(P->P.predStrings.isEmpty())
			)
			{
			if(this.removeVariantIfNoMoreAttribute) continue;
			variant_has_one_matching_pred=false;
			}
		
		/* all handlers failed */
		if(variant_has_one_matching_pred && 
			!this.removeUnusedAttribute &&
			 predictionHandlers.stream().allMatch(P->!P.keepFlag))
			{
			variant_has_one_matching_pred=false;
			}
		
		
		if(!StringUtil.isBlank(this.filterIn))
			{
			if(variant_has_one_matching_pred){
				vcb.filter(this.filterIn);
				}
			else if( !ctx.filtersWereApplied()) {
				vcb.passFilters();
				}
			}
		else  if(!StringUtil.isBlank(this.filterOut)) {
			if(variant_has_one_matching_pred && !ctx.filtersWereApplied()) {
				vcb.passFilters();
				}
			else if(!variant_has_one_matching_pred) {
				vcb.filter(this.filterOut);
				}
			}
		else if(!variant_has_one_matching_pred)
			{
			continue;
			}
		
		final VariantContext ctx3;
		if(GT_FILTER_RESET_TO_NOCALL.equals(this.filterGenotypesStr))
		 	{
			ctx3 =this.recalculator.apply(vcb.make());
		 	}
		 else
		 	{
			ctx3= vcb.make(); 
		 	}
		return ctx3;
		}

				
			
	private boolean hasUserTemLabel(final Collection<SequenceOntologyTree.Term> ctxTerms)
		{
		if(ctxTerms==null || ctxTerms.isEmpty()) return false;
		return ctxTerms.
				stream().
				anyMatch(S->this.user_terms.contains(S));
		}				
				
	private boolean isRecodingGenotypes() {
		return !StringUtil.isBlank(this.filterGenotypesStr);
	}
	
	
	public void compileTerms(final Collection<String> userTermsAsString) {			
			/* map term as string to Term instances */
			final Set<SequenceOntologyTree.Term> tmpSet1 = new HashSet<>();
		    userTermsAsString.
		    	stream().
				map(S->S.trim()).
				filter(L->!StringUtils.isBlank(L)).
				flatMap(S->Arrays.stream(S.split("[,;& \t]"))).
				filter(L->!StringUtils.isBlank(L)).
				map(ACN->{
					SequenceOntologyTree.Term T = getSequenceOntologyTree().getTermByAcn(ACN);
					if(T==null)
						{
						T = getSequenceOntologyTree().getTermByLabel(ACN);
						}
					if(T==null)
						{
						throw new JvarkitException.UserError("Unknown SO:Accession/label \""+ACN+"\"");
						}
					return T;					
					}).
				forEach(t->{
					tmpSet1.add(t);
					/* reasoning enabled */
					if(!this.disableReasoning) tmpSet1.addAll(t.getAllDescendants());					
				});

		    /* inverse logic if needed */
			if(this.invertSoTerms)
				{
				final Set<SequenceOntologyTree.Term> tmpSet2 = new HashSet<>(getSequenceOntologyTree().getTerms());
				tmpSet2.removeAll(tmpSet1);
				this.user_terms.addAll(tmpSet2);
				}
			else
				{
				this.user_terms.addAll(tmpSet1);
				}
			}
	
	@Override
	public void fillHeader(VCFHeader header) {
		final String termlist = String.join(", ",this.user_terms.stream().
				map(S->S.getAcn()+"("+S.getLabel()+")").
				collect(Collectors.toSet()))
				;
		
		if(!StringUtil.isBlank(this.filterIn) && !StringUtil.isBlank(this.filterOut)) {
			throw new IllegalArgumentException("Option filterIn && filterOut both defined.");
			}

		
		if(!StringUtils.isBlank(this.filterIn)) {
			header.addMetaDataLine(new VCFFilterHeaderLine(this.filterIn,
					"Variant having SO terms:"+ termlist));
			}
		if(!StringUtils.isBlank(this.filterOut)) {
			header.addMetaDataLine(new VCFFilterHeaderLine(this.filterOut,
					"Variant non having SO terms :" + termlist));
			}
	
		if(!StringUtil.isBlank(this.filterGenotypesStr) &&
			!GT_FILTER_RESET_TO_NOCALL.equals(this.filterGenotypesStr))
			{
			header.addMetaDataLine(new VCFFormatHeaderLine(
					VCFConstants.GENOTYPE_FILTER_KEY,
					1,VCFHeaderLineType.String,
					"Genotype was filterered by vcffilterso : "+this.filterGenotypesStr
					));
			}
	
		AbstractPredictionHandler ph = new VepPredictionHandler(header);
		if(ph.isValid()) predictionHandlers.add(ph);
		ph = new SnpEffPredictionHandler(header);
		if(ph.isValid()) predictionHandlers.add(ph);
		ph = new AnnPredictionHandler(header);
		if(ph.isValid()) predictionHandlers.add(ph);
		ph = new BcftoolsCsqPredictionHandler(header);
		if(ph.isValid()) predictionHandlers.add(ph);
		}
	
	
	@Override
	public List<VariantContext> annotate(VariantContext ctx) throws IOException {
		return Collections.singletonList(annotate1(ctx));
		}
	
	@Override
	public void close() {
		
		}
	}
