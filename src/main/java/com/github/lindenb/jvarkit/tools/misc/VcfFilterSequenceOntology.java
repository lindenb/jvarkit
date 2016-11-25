/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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


History:
* 2014 creation
* 2015 knime interface

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.MyPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;


/**
 * VCF filter on Sequence Ontology
 * @author lindenb
 *
 */
public class VcfFilterSequenceOntology
	extends AbstractVcfFilterSequenceOntology
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfFilterSequenceOntology.class);
	
	
	/* all sequence terms */
	private final Set<SequenceOntologyTree.Term> user_terms=new HashSet<SequenceOntologyTree.Term>();
	/* SO Treem */
	private  SequenceOntologyTree sequenceOntologyTree=SequenceOntologyTree.getInstance();
	
	
	/* public : knime needs this*/
	public VcfFilterSequenceOntology()
		{
		}
	
	
	private void addVariant(final VariantContextWriter w,final VariantContext ctx,boolean keep)
		{
		if(super.isInvert()) keep=!keep;
		if(!this.filterIn.isEmpty())
			{
			if(keep){
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.filter(this.filterIn);
				w.add(vcb.make());
				}
			else if( !ctx.filtersWereApplied()) {
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.passFilters();
				w.add(vcb.make());
				}
			else
				{
				w.add(ctx);
				}
			}
		else  if(!this.filterOut.isEmpty()) {
			if(keep && !ctx.filtersWereApplied()) {
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.passFilters();
				w.add(vcb.make());
				}
			else if(keep){
				w.add(ctx);
				}
			else
				{
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.filter(this.filterOut);
				w.add(vcb.make());
				}
			}
		else
			{
			if(keep) {
				w.add(ctx);
				} else
				{
					/* don't print */
				}
			}
		}

	
	private boolean hasUserTem(final Set<SequenceOntologyTree.Term> ctxTerms)
		{
		return !Collections.disjoint(ctxTerms, this.user_terms);
		}
	
	@Override
	/* public for knime */ public Collection<Throwable> doVcfToVcf(final String inputName,final  VcfIterator in,final  VariantContextWriter out) {
		try {
			final VCFHeader header=in.getHeader();
			final VCFHeader h2=new VCFHeader(header);
			addMetaData(h2);
			
			out.writeHeader(h2);

			final VepPredictionParser vepParser=new VepPredictionParser(header).sequenceOntologyTree(this.sequenceOntologyTree);
			final SnpEffPredictionParser snpEffparser= new SnpEffPredictionParser(header).sequenceOntologyTree(this.sequenceOntologyTree);
			final MyPredictionParser myPredParser= new MyPredictionParser(header).sequenceOntologyTree(this.sequenceOntologyTree);
			final AnnPredictionParser annPredParser= new AnnPredictionParser(header).sequenceOntologyTree(this.sequenceOntologyTree);
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
			
			while(in.hasNext() )
				{	
				final  VariantContext ctx=progress.watch(in.next());
				boolean keep=false;
				List<Object> snpEffPredStrings = null;
				List<Object> vepPredStrings = null;
				List<Object> myPredStrings = null;
				List<Object> annPredStrings = null;
				
				/* handle SNP EFF */
				if(ctx.hasAttribute(snpEffparser.getTag())) {
					snpEffPredStrings = new ArrayList<>(ctx.getAttributeAsList(snpEffparser.getTag()));
					int x=0;
					while(x<snpEffPredStrings.size()){
						final Object predStr = snpEffPredStrings.get(x);
						final SnpEffPredictionParser.SnpEffPrediction pred = snpEffparser.parseOnePrediction(predStr);
						if(pred!=null && hasUserTem(pred.getSOTerms()))
							{
							keep=true;
							++x;
							}
						else
							{
							snpEffPredStrings.remove(x);
							}
						}
					}
				
				/* handle VEP */
				if(ctx.hasAttribute(vepParser.getTag()) && (super.removeUnusedAttribute || !keep)) {
					vepPredStrings = new ArrayList<>(ctx.getAttributeAsList(vepParser.getTag()));
					int x=0;
					while(x<vepPredStrings.size()){
						final Object predStr = vepPredStrings.get(x);
						final VepPredictionParser.VepPrediction pred = vepParser.parseOnePrediction(ctx, predStr);
						if(pred!=null && hasUserTem(pred.getSOTerms()))
							{
							keep=true;
							++x;
							}
						else
							{
							vepPredStrings.remove(x);
							}
						}
					}
				/* handle MyPredictionParser */
				if(ctx.hasAttribute(myPredParser.getTag()) && (super.removeUnusedAttribute || !keep)) {
					myPredStrings = new ArrayList<>(ctx.getAttributeAsList(myPredParser.getTag()));
					int x=0;
					while(x<myPredStrings.size()){
						final Object predStr = myPredStrings.get(x);
						final MyPredictionParser.MyPrediction pred = myPredParser.parseOnePrediction(predStr);
						if(pred!=null && hasUserTem(pred.getSOTerms()))
							{
							keep=true;
							++x;
							}
						else
							{
							myPredStrings.remove(x);
							}
						}
					}
				
				/* handle ANN */
				if(ctx.hasAttribute(annPredParser.getTag()) && (super.removeUnusedAttribute || !keep)) {
					annPredStrings = new ArrayList<>(ctx.getAttributeAsList(annPredParser.getTag()));
					int x=0;
					while(x<annPredStrings.size()){
						final Object predStr = annPredStrings.get(x);
						final AnnPredictionParser.AnnPrediction pred = annPredParser.parseOnePrediction(predStr);
						if(pred!=null && hasUserTem(pred.getSOTerms()))
							{
							keep=true;
							++x;
							}
						else
							{
							annPredStrings.remove(x);
							}
						}
					}
				
				if(super.removeUnusedAttribute)
					{
					final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
					int count=0;
					if(snpEffPredStrings!=null) {
						vcb.rmAttribute(snpEffparser.getTag());
						vcb.attribute(snpEffparser.getTag(), snpEffPredStrings);
						count += snpEffPredStrings.size();
						}
					if(vepPredStrings!=null) {
						vcb.rmAttribute(vepParser.getTag());
						vcb.attribute(vepParser.getTag(), vepPredStrings);
						count += vepPredStrings.size();
						}
					if(myPredStrings!=null) {
						vcb.rmAttribute(myPredParser.getTag());
						vcb.attribute(myPredParser.getTag(), myPredStrings);
						count += myPredStrings.size();
						}
					if(annPredStrings!=null) {
						vcb.rmAttribute(annPredParser.getTag());
						vcb.attribute(annPredParser.getTag(), annPredStrings);
						count += annPredStrings.size();
						}
					
					if( count == 0  && this.removeIfNoMoreAttribute) 
						{
						addVariant(out,vcb.make(),false);
						}
					else
						{
						addVariant(out,vcb.make(),true);
						}
					}
				else
					{
					addVariant(out,ctx,keep);
					}
				}
			progress.finish();
			return RETURN_OK;
			}
		finally
			{
			}
		}
	
	
	
	private void parseAccessionsFile(final File f) throws IOException
		{
		final BufferedReader in=IOUtils.openFileForBufferedReading(f);
		String line;
		while((line=in.readLine())!=null)
			{
			if(line.startsWith("#")) continue;
			line=line.trim();
			if(line.trim().isEmpty()) continue;
			super.userTermsAsString.add(line);
			}
		in.close();
		}
	
	@Override
	protected VCFHeader addMetaData(final VCFHeader header) {
		final String termlist = String.join(", ",this.user_terms.stream().map(S->S.getAcn()).collect(Collectors.toSet()));
		if(!super.filterIn.isEmpty()) {
			header.addMetaDataLine(new VCFFilterHeaderLine(super.filterIn,
					"Variant having SO terms:"+ termlist));
			}
		if(!super.filterOut.isEmpty()) {
			header.addMetaDataLine(new VCFFilterHeaderLine(super.filterOut,
					"Variant non having SO terms :" + termlist));
			}
		return super.addMetaData(header);
		}

	
	@Override
	public Collection<Throwable> initializeKnime() {
		if(this.removeIfNoMoreAttribute) {
			 super.removeUnusedAttribute = true;
		}
		
		if(super.invert && super.removeUnusedAttribute) {
			return wrapException("Option -"+OPTION_INVERT+" cannot be used when Option -"+OPTION_REMOVEUNUSEDATTRIBUTE);
		}

		
		
		if(!super.filterIn.isEmpty() && !super.filterOut.isEmpty()) {
			return wrapException("Option -"+OPTION_FILTERIN+"  and  -"+OPTION_FILTEROUT+" both defined.");
		}
		if((super.invert) && (!super.filterIn.isEmpty() || !super.filterOut.isEmpty())) {
			return wrapException("Option -"+OPTION_INVERT+"cannot be used when Option -"+OPTION_FILTERIN+" or  -"+OPTION_FILTEROUT+" is defined.");
		}
		
		
		if( !(super.owluri==null || super.owluri.trim().isEmpty()) ) {
			LOG.info("loading so tree from "+super.owluri);
			try
				{
				this.sequenceOntologyTree = SequenceOntologyTree.fromUri(super.owluri.trim());
				}
			catch (final IOException e)
				{
				return wrapException(e);
				}
			LOG.info("Done loading SO.");
			}
		
		final boolean reasoning = !super.disableReasoning;
		if(super.userAcnFile!=null)
			{
			try {
				this.parseAccessionsFile(super.userAcnFile);
			} catch (Exception e) {
				return wrapException(e);
				}
			}
		
		
		for(String acn: super.userTermsAsString)
			{
			acn=acn.trim();
			if(acn.isEmpty()) continue;
			final SequenceOntologyTree.Term t= this.sequenceOntologyTree.getTermByAcn(acn);
			if(t==null)
				{
				return wrapException("Unknown SO:Accession \""+acn+"\"");
				}
			this.user_terms.add(t);
			if(reasoning) this.user_terms.addAll(t.getAllDescendants());
			}
		if(this.user_terms.isEmpty())
			{
			LOG.warn("No SO: term found ");
			}
		LOG.info("Will be using :"+this.user_terms.toString());
		
		this.userTermsAsString.clear();//we don't need this anymore
		return super.initializeKnime();
		}
	@Override
	public void disposeKnime() {
		super.userTermsAsString.clear();
		this.user_terms.clear();
		super.disposeKnime();
		}
		
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		try {
			if(super.showList)
				{
				PrintWriter pw=super.openFileOrStdoutAsPrintWriter();
				for(SequenceOntologyTree.Term t:sequenceOntologyTree.getTerms())
					{
					pw.println(t.getAcn()+"\t"+t.getLabel());
					}
				pw.close();
				return RETURN_OK;
				}
			return doVcfToVcf(inputName);
		} catch (Exception e) {
			}
		finally {
			
		}
		
		return null;
		}
	
	public static void main(String[] args)
		{
		new VcfFilterSequenceOntology().instanceMainWithExit(args);
		}
	}
