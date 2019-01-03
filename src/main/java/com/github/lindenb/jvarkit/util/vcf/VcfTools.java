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
package com.github.lindenb.jvarkit.util.vcf;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.github.lindenb.jvarkit.annotproc.IncludeSourceInJar;
import com.github.lindenb.jvarkit.tools.vcftrios.DeNovoDetector;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/** bean to be injected in various contexts like javascript context, etc... */
@IncludeSourceInJar
public class VcfTools {
public static final Logger LOG = Logger.build(VcfTools.class).make();

private VCFHeader header=null;
private SnpEffPredictionParser snpEffPredictionParser=null;
private VepPredictionParser vepPredictionParser=null;
private AnnPredictionParser annPredictionParser=null;
private final DeNovoDetector deNovoDetector = new DeNovoDetector();
public VcfTools() {
	init(null);
	}
public VcfTools(final VCFHeader header) {
	init(header);
	}

public VCFHeader getHeader() {
	return header;
}

public void init(final VCFHeader header) {
	this.header=header;
	this.snpEffPredictionParser=new SnpEffPredictionParserFactory().header(header).get();
	this.vepPredictionParser=new VepPredictionParserFactory().header(header).get();
	this.annPredictionParser=new AnnPredictionParserFactory(header).get();	
	}

public VepPredictionParser getVepPredictionParser() {
	return vepPredictionParser;
	}

public SnpEffPredictionParser getSnpEffPredictionParser() {
	return snpEffPredictionParser;
	}

public AnnPredictionParser getAnnPredictionParser() {
	return annPredictionParser;
	}
public void initSnpEffParser(final String definition)
	{
	failIf(definition==null || definition.trim().isEmpty(),"SnpEff definition is empty");
	final VCFHeader header=new VCFHeader();
	final VCFInfoHeaderLine info=new VCFInfoHeaderLine(
			VepPredictionParser.getDefaultTag(),
			VCFHeaderLineCount.UNBOUNDED,
			VCFHeaderLineType.String,
			definition
			);
	header.addMetaDataLine(info);
	this.vepPredictionParser=new VepPredictionParserFactory(header).get();
	}
protected static void failIf(boolean testFailed, String msg) {
	if (!testFailed)
		return;
	LOG.error("FAILURE: " + msg);
	throw new RuntimeException("The following test failed: " + msg);
	}
public SequenceOntologyTree getSequenceOntologyTree() {
	return  SequenceOntologyTree.getInstance();
}

public List<AnnPredictionParser.AnnPrediction> getAnnPredictions(final VariantContext ctx) {
	if(this.annPredictionParser==null) return Collections.emptyList();
	return this.annPredictionParser.getPredictions(ctx);
	}

public List<VepPredictionParser.VepPrediction> getVepPredictions(final VariantContext ctx) {
	if(this.getVepPredictionParser()==null) return Collections.emptyList();
	return this.getVepPredictionParser().getPredictions(ctx);
	}

public List<SnpEffPredictionParser.SnpEffPrediction> getSnpEffPredictions(final VariantContext ctx) {
	if(this.getSnpEffPredictionParser()==null) return Collections.emptyList();
	return this.getSnpEffPredictionParser().getPredictions(ctx);
	}


/** return true if variant has any prediction with a SO term (or its children) with this label */
public boolean hasSequenceOntologyLabel(final VariantContext ctx,final String lbl)
	{
	if(lbl==null) return false;
	final SequenceOntologyTree.Term t= this.getSequenceOntologyTree().getTermByLabel(lbl);
	if(t==null) LOG.warning("don't know SO.label "+lbl);
	return hasSequenceOntologyTerm(ctx,t);
	}
/** return true if variant has any prediction with a SO term (or its children) with this accession */
public boolean hasSequenceOntologyAccession(final VariantContext ctx,final String acn)
	{
	if(acn==null) return false;
	final SequenceOntologyTree.Term t= this.getSequenceOntologyTree().getTermByAcn(acn);
	if(t==null) LOG.warning("don't know SO.acn "+acn);
	return hasSequenceOntologyTerm(ctx,t);
	}

/** return true if variant has any prediction with a SO term (or its children) */
public boolean hasSequenceOntologyTerm(final VariantContext ctx,final SequenceOntologyTree.Term t)
	{
	if(t==null) return false;
	final Set<SequenceOntologyTree.Term> children=t.getAllDescendants();
	for(final AnnPredictionParser.AnnPrediction a: getAnnPredictions(ctx)) {
		if(!Collections.disjoint(a.getSOTerms(),children)) return true;
		}
	for(final VepPredictionParser.VepPrediction a: getVepPredictions(ctx)) {
		if(!Collections.disjoint(a.getSOTerms(),children)) return true;
		}
	for(final SnpEffPredictionParser.SnpEffPrediction a: getSnpEffPredictions(ctx)) {
		if(!Collections.disjoint(a.getSOTerms(),children)) return true;
		}
	
	
	return false;
	}

public boolean isMendelianIncompatibility(final Genotype child,final Genotype parent)
	{
	if(child==null || parent==null) return false;
	return getDeNovoDetector().test(null, parent, child)!=null;
	}

public DeNovoDetector getDeNovoDetector() {
	return this.deNovoDetector;
}

public boolean isMendelianIncompatibility(final Genotype child,final Genotype father,final Genotype mother)
	{
	return getDeNovoDetector().test(father, mother, child)!=null;
	}
/** return true if there is a mendelian problem with the children */
public boolean isMendelianIncompatibility(final VariantContext ctx,final Pedigree.Person child)
	{
	final Genotype gc= ctx.getGenotype(child.getId());
	if(gc==null) return false;
	final Genotype gf= child.hasFather()?ctx.getGenotype(child.getFatherId()):null;
	final Genotype gm= child.hasMother()?ctx.getGenotype(child.getMotherId()):null;
	if(gf==null && gm == null) return false;
	return getDeNovoDetector().test(ctx,gf,gm,gc)!=null;
	}

}
