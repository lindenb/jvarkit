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
package com.github.lindenb.jvarkit.tools.vcfsparql;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.BiPredicate;
import java.util.function.Predicate;

import org.apache.jena.datatypes.xsd.XSDDatatype;
import org.apache.jena.graph.Node;
import org.apache.jena.graph.NodeFactory;
import org.apache.jena.graph.Triple;
import org.apache.jena.graph.TripleIterator;
import org.apache.jena.graph.impl.GraphBase;
import org.apache.jena.util.iterator.ExtendedIterator;
import org.apache.jena.util.iterator.NiceIterator;
import org.apache.jena.vocabulary.RDF;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser.AnnPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;


public class VariantGraph extends GraphBase {
	static final String NS="jvarkit:";
	private final Node RDF_TYPE = RDF.type.asNode();
	private final Path variantFile;
	private final Map<String,Node> cache = new HashMap<>();
	private Predicate<VariantContext> acceptVariant = V->true;
	private BiPredicate<VariantContext,Genotype> acceptGenotype = (V,G)->true;
	private Locatable interval = null;
	private AnnPredictionParser annParser = null;
	
	public VariantGraph(final Path variantFile) {
		this.variantFile = variantFile;
		try(final VCFFileReader r=new VCFFileReader(variantFile,false)) {
			final VCFHeader header  = r.getFileHeader();
			final AnnPredictionParserFactory f = new AnnPredictionParserFactory(header);
			this.annParser = f.get();
			}
		}
	@Override
	protected ExtendedIterator<Triple> graphBaseFind(final Triple triple) {
		final TripleVcfIterator iter0 = new TripleVcfIterator();
		final TripleAdaptorIterator iter1 = new TripleAdaptorIterator(iter0);
		if(triple==null) return iter1;
		return iter1.filterKeep((T)->triple.matches(T));
		}
	
	public void setVariantFilter(final Predicate<VariantContext> pred) {
		this.acceptVariant = pred==null?V->true:pred;
	}
	public void setGenotypeFilter(final BiPredicate<VariantContext,Genotype> pred) {
		this.acceptGenotype = pred==null?(V,G)->true:pred;
	}
	
	public void setInterval(Locatable interval) {
		this.interval = interval;
		}
	
	@Override
	public void close() {
		super.close();
		}
	
	
	private Node mkNode(final String local) {
		final String NS="jvarkit:";
		return NodeFactory.createURI(NS+local);
		}
	
	private Node cachedNode(final String local) {
		Node n = this.cache.get(local);
		if(n==null) {
			n = mkNode(local);
			this.cache.put(local, n);
			}
		return n;
		}
	

	
	protected List<Triple> variantToTriple(final long nLine,final VariantContext ctx)
		{
		final List<Triple> triples = new ArrayList<>();
		final Node ctxId = NodeFactory.createBlankNode(String.format("var%06d", nLine));
		
		triples.add(Triple.create(ctxId, RDF_TYPE,cachedNode("Variant")));
		triples.add(Triple.create(ctxId, cachedNode("file"),NodeFactory.createLiteral(this.variantFile.toString())));
		triples.add(Triple.create(ctxId, RDF_TYPE,cachedNode(ctx.getType().name())));
		triples.add(Triple.create(ctxId,cachedNode("contig"),NodeFactory.createLiteral(ctx.getContig())));
		triples.add(Triple.create(ctxId,cachedNode("start"),NodeFactory.createLiteral(String.valueOf(ctx.getStart()),XSDDatatype.XSDint)));
		triples.add(Triple.create(ctxId,cachedNode("end"),NodeFactory.createLiteral(String.valueOf(ctx.getEnd()),XSDDatatype.XSDint)));
		if(ctx.hasID()) {
			triples.add(Triple.create(ctxId,cachedNode("id"),NodeFactory.createLiteral(ctx.getID())));
		}
		//REF + ALT
		for(final Allele allele : ctx.getAlleles()) {
			triples.add(Triple.create(ctxId,cachedNode(allele.isReference()?"ref":"alt"),
				NodeFactory.createLiteral(allele.getDisplayString())));
			}
		// QUAL
		if(ctx.hasLog10PError())
			{
			
			}
		for(final String filter : ctx.getFilters()) {
			triples.add(Triple.create(ctxId,cachedNode("filter"),
				NodeFactory.createLiteral(filter)));
			}
		
		if(this.annParser!=null && this.annParser.isValid()) {
			for(final AnnPrediction p:this.annParser.getPredictions(ctx)) {
				final Node predIdx = NodeFactory.createBlankNode();
				triples.add(Triple.create(predIdx, RDF_TYPE,cachedNode("SnpEff")));
				triples.add(Triple.create(ctxId, cachedNode("prediction"),predIdx));
				if(!StringUtils.isBlank(p.getFeatureId())) {
					triples.add(Triple.create(predIdx, cachedNode("feature-id"),NodeFactory.createLiteral(p.getFeatureId())));
					}
				if(p.getPutativeImpact()!=null) {
					triples.add(Triple.create(predIdx, cachedNode("impact"),NodeFactory.createLiteral(p.getPutativeImpact().name())));
					}
				if(!StringUtils.isBlank(p.getGeneName())) {
					triples.add(Triple.create(predIdx, cachedNode("gene-name"),NodeFactory.createLiteral(p.getGeneName())));
					}
				for(final String so:p.getSOTermsStrings())
					{
					triples.add(Triple.create(predIdx, cachedNode("so-term"),NodeFactory.createLiteral(so)));
					}
				}
		}
		
		for(int i=0;i< ctx.getNSamples();++i)
			{
			final Genotype gt=ctx.getGenotype(i);
			if(!this.acceptGenotype.test(ctx, gt)) continue;
			
			final Node gtId = NodeFactory.createBlankNode(String.format("var%06d.gt%03d", nLine,(i+1)));
			
			triples.add(Triple.create(gtId, RDF_TYPE,cachedNode("Genotype")));
			triples.add(Triple.create(ctxId, cachedNode("genotype"),gtId));
			triples.add(Triple.create(ctxId, cachedNode("gtype"),NodeFactory.createLiteral(gt.getType().name())));
			triples.add(Triple.create(gtId, cachedNode("sample"),NodeFactory.createLiteral(gt.getSampleName())));
			if(gt.hasDP())
				{
				triples.add(Triple.create(gtId,cachedNode("dp"),NodeFactory.createLiteral(String.valueOf(gt.getDP()),XSDDatatype.XSDint)));
				}
			if(gt.hasGQ())
				{
				triples.add(Triple.create(gtId,cachedNode("gq"),NodeFactory.createLiteral(String.valueOf(gt.getDP()),XSDDatatype.XSDint)));
				}
			for(final Allele allele : gt.getAlleles()) {
				triples.add(Triple.create(gtId,cachedNode("allele"),
					NodeFactory.createLiteral(allele.getDisplayString())));
				}
			}
		
		
		return triples;
		}
	
	private class TripleAdaptorIterator
	extends NiceIterator<Triple>
	implements TripleIterator
		{
		final TripleVcfIterator delegate;
		TripleAdaptorIterator(final TripleVcfIterator delegate) {
			this.delegate = delegate;
			}
		@Override
		public Triple nextTriple() {
			return this.next();
			}
		@Override
		public boolean hasNext() {
			return this.delegate.hasNext();
			}
		@Override
		public Triple next() {
			return this.delegate.next();
			}
		}
		
	private class TripleVcfIterator
		extends AbstractIterator<Triple>
		implements CloseableIterator<Triple>
		{
		private final VCFFileReader vcfFileReader;
		private final CloseableIterator<VariantContext> vcfIterator;
		private final List<Triple> buffer=new ArrayList<>();
		private long nLines=0;
		TripleVcfIterator() {
			this.vcfFileReader = new VCFFileReader(VariantGraph.this.variantFile,
					interval==null?false:true);
			this.vcfIterator = 
					interval==null?
					this.vcfFileReader.iterator():
					this.vcfFileReader.query(interval);
			}
		
		@Override
		protected Triple advance() {
			while(this.buffer.isEmpty())
				{
				if(!this.vcfIterator.hasNext()) return null;
				final VariantContext ctx = this.vcfIterator.next();
				++nLines;
				if(!acceptVariant.test(ctx)) continue;
				if(!this.vcfIterator.hasNext()) this.close();
				this.buffer.addAll(variantToTriple(nLines,ctx));
				}
			return this.buffer.remove(this.buffer.size()-1);
			}
		
		@Override
		public void close() {
			this.buffer.clear();
			this.vcfIterator.close();
			this.vcfFileReader.close();
			}
		}
	
	}
