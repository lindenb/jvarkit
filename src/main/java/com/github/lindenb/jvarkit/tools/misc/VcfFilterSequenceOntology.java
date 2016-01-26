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
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
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
	private Set<SequenceOntologyTree.Term> user_terms=new HashSet<SequenceOntologyTree.Term>();
	/* inverse result */
	private final SequenceOntologyTree sequenceOntologyTree=SequenceOntologyTree.getInstance();
	
	
	/* public : knime needs this*/
	public VcfFilterSequenceOntology()
		{
		}
	
	
	private boolean hasUserTem(final Set<SequenceOntologyTree.Term> ctxTerms)
		{
		for(final SequenceOntologyTree.Term ctxTerm:ctxTerms)
			{
			if(this.user_terms.contains(ctxTerm))
				{
				return true;
				}
			}
		return false;
		}
	
	@Override
	/* public for knime */ public Collection<Throwable> doVcfToVcf(final String inputName,final  VcfIterator in,final  VariantContextWriter out) {
		try {
			VCFHeader header=in.getHeader();
			VCFHeader h2=new VCFHeader(header);
			addMetaData(h2);
			
			out.writeHeader(h2);

			final VepPredictionParser vepParser=new VepPredictionParser(header);
			final SnpEffPredictionParser snpEffparser=new SnpEffPredictionParser(header);
			final MyPredictionParser myPredParser=new MyPredictionParser(header);
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			while(in.hasNext() )
				{	
				VariantContext ctx=progress.watch(in.next());
				boolean keep=false;
								
				for(SnpEffPredictionParser.SnpEffPrediction pred:snpEffparser.getPredictions(ctx))
					{
					if(hasUserTem(pred.getSOTerms())) { keep=true; break;}
					}
				if(!keep)
					{
					for(VepPredictionParser.VepPrediction pred:vepParser.getPredictions(ctx))
						{
						if(hasUserTem(pred.getSOTerms())) { keep=true; break;}
						}
					}
				if(!keep)
					{
					for(MyPredictionParser.MyPrediction pred:myPredParser.getPredictions(ctx))
						{
						if(hasUserTem(pred.getSOTerms())) { keep=true; break;}
						}
					}
				if(super.isInvert() ) keep=!keep;
				if(keep)
					{
					out.add(ctx);
					}
				if(out.checkError()) break;
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
		BufferedReader in=IOUtils.openFileForBufferedReading(f);
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
	public Collection<Throwable> initializeKnime() {
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
			SequenceOntologyTree.Term t=sequenceOntologyTree.getTermByAcn(acn);
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
		
		this.userTermsAsString.clear();;//we don't need this anymore
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
