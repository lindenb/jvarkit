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

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import com.github.lindenb.jvarkit.util.command.Command;
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
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(VcfFilterSequenceOntology.class);

	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractVcfFilterSequenceOntology.AbstractVcfFilterSequenceOntologyCommand
		{		
	/* all sequence terms */
		private Set<SequenceOntologyTree.Term> user_terms=new HashSet<SequenceOntologyTree.Term>();
		/* inverse result */
		private final SequenceOntologyTree sequenceOntologyTree=SequenceOntologyTree.getInstance();
	
		private boolean hasUserTem(Set<SequenceOntologyTree.Term> ctxTerms)
			{
			for(SequenceOntologyTree.Term ctxTerm:ctxTerms)
				{
				if(this.user_terms.contains(ctxTerm))
					{
					return true;
					}
				}
			return false;
			}
		@Override
		protected Collection<Throwable> doVcfToVcf(String inputName,
				VcfIterator in, VariantContextWriter out) throws IOException {
		try {
			VCFHeader header=in.getHeader();
			VCFHeader h2=new VCFHeader(header);
			addMetaData(h2);
			out.writeHeader(h2);

			final VepPredictionParser vepParser=new VepPredictionParser(header);
			final SnpEffPredictionParser snpEffparser=new SnpEffPredictionParser(header);
			final MyPredictionParser myPredParser=new MyPredictionParser(header);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
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
				if(super.inverse_result ) keep=!keep;
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
	
	
	@Override
	public Collection<Throwable> initializeKnime()
		{
		if(super.termFile!=null)
			{
			try {
				super.userTermsAsString.addAll(super.readStringSet(super.termFile));
			} catch (Exception e) {
				return wrapException(e);
				}
			}
		
		
		for(String acn: this.userTermsAsString)
			{
			acn=acn.trim();
			if(acn.isEmpty()) continue;
			SequenceOntologyTree.Term t=sequenceOntologyTree.getTermByAcn(acn);
			if(t==null)
				{
				return wrapException("Unknown SO:Accession \""+acn+"\"");
				}
			this.user_terms.add(t);
			if(!disable_reasoning) this.user_terms.addAll(t.getAllDescendants());
			}
		if(this.user_terms.isEmpty())
			{
			LOG.warn("No SO: term found ");
			}
		LOG.info("Will be using :"+this.user_terms.toString());
		
		this.userTermsAsString=null;//we don't need this anymore
		return super.initializeKnime();
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		if(super.just_list)
			{
			PrintWriter out=null;
			try {
				out= openFileOrStdoutAsPrintWriter();
				for(SequenceOntologyTree.Term t:sequenceOntologyTree.getTerms())
					{
					out.println(t.getAcn()+"\t"+t.getLabel());
					}
			} catch (Exception e) {
				return wrapException(e);
				}
			finally
				{
				CloserUtil.close(out);
				}
			}
		return doVcfToVcf(inputName);
		}
	}
	

	public static void main(String[] args)
		{
		new VcfFilterSequenceOntology().instanceMainWithExit(args);
		}
	}
