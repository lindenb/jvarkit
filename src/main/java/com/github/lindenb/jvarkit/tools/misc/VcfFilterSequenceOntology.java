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
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.knime.KnimeApplication;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
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
	extends AbstractCommandLineProgram
	implements KnimeApplication
	{
	private Set<String> userTermsAsString=new HashSet<String>();
	/* enable reasoning */
	private boolean reasoning=true;
	/* all sequence terms */
	private Set<SequenceOntologyTree.Term> user_terms=new HashSet<SequenceOntologyTree.Term>();
	/* inverse result */
	private boolean inverse_result=false;
	private final SequenceOntologyTree sequenceOntologyTree=SequenceOntologyTree.getInstance();
	/** output file : stdout if null */
	private File outputFile=null;
	/** number of variants filterer */
	private int countFilteredVariants=0;
	
	
	/* public : knime needs this*/
	public VcfFilterSequenceOntology()
		{
		}
	

	@Override
	public String getProgramDescription() {
		return "Filter a VCF file annotated with SNPEff or VEP with terms from Sequence-Ontology." +
				" Reasoning : Children of user's SO-terms will be also used.";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfFilterSequenceOntology";
		}
	
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
	public int executeKnime(List<String> args)
		{
		VcfIterator vcfIn=null;
		try
			{
			if(args.isEmpty())
				{
				vcfIn = VCFUtils.createVcfIteratorStdin();
				}
			else if(args.size()==1)
				{
				vcfIn= VCFUtils.createVcfIterator(args.get(0));
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			this.filterVcfIterator(vcfIn);
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(vcfIn);
			}
		
		}
	
	private void filterVcfIterator(VcfIterator in) throws IOException
		{
		this.countFilteredVariants=0;
		VariantContextWriter out = null;
		try {
			VCFHeader header=in.getHeader();
			VCFHeader h2=new VCFHeader(header);
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));

			
			if(getOutputFile()==null)
				{
				out = VCFUtils.createVariantContextWriterToStdout();
				}
			else
				{
				info("opening vcf writer to "+getOutputFile());
				out = VCFUtils.createVariantContextWriter(getOutputFile());
				}
			out.writeHeader(h2);

			final VepPredictionParser vepParser=new VepPredictionParser(header);
			final SnpEffPredictionParser snpEffparser=new SnpEffPredictionParser(header);
			final MyPredictionParser myPredParser=new MyPredictionParser(header);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			while(in.hasNext() )
				{	
				this.checkKnimeCancelled();
				VariantContext ctx=in.next();
				boolean keep=false;
				
				
				progress.watch(ctx.getChr(), ctx.getStart());
				
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
				if(isInverseResult() ) keep=!keep;
				if(keep)
					{
					this.countFilteredVariants++;
					out.add(ctx);
					}
				}
			progress.finish();
			}
		finally
			{
			CloserUtil.close(out);
			out=null;
			}
		}
	
	/** return the number of variant kept during last invocation of this program */
	public int geVariantCount()
		{
		return countFilteredVariants;
		}
	
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -A (SO:ACN). add this SO:ACN");
		out.println(" -f (filename). Tab delimited file of SO accession numbers");
		out.println(" -S list the available SO accession and exit.");
		out.println(" -v invert selection.");
		out.println(" -d disable reasoning, don't use term's children.");
		out.println(" -o (output-file) out (default: stdout).");
		super.printOptions(out);
		}
	
	private void parseAccessionsFile(File f) throws IOException
		{
		BufferedReader in=IOUtils.openFileForBufferedReading(f);
		String line;
		while((line=in.readLine())!=null)
			{
			if(line.startsWith("#")) continue;
			line=line.trim();
			if(line.trim().isEmpty()) continue;
			addTerm(line);
			}
		in.close();
		}
	
	@Override
	public void checkKnimeCancelled() {
		
		}
	
	@Override
	public void disposeKnime() {
		
		}
	@Override
	public int initializeKnime()
		{
		for(String acn: this.userTermsAsString)
			{
			acn=acn.trim();
			if(acn.isEmpty()) continue;
			SequenceOntologyTree.Term t=sequenceOntologyTree.getTermByAcn(acn);
			if(t==null)
				{
				error("Unknown SO:Accession \""+acn+"\"");
				return -1;
				}
			this.user_terms.add(t);
			if(reasoning) this.user_terms.addAll(t.getAllDescendants());
			}
		if(this.user_terms.isEmpty())
			{
			warning("No SO: term found ");
			}
		info("Will be using :"+this.user_terms.toString());
		
		this.userTermsAsString=null;//we don't need this anymore
		return 0;
		}
	
	public void addTerm(String acn)
		{
		this.userTermsAsString.add(acn);
		}
	
	public void setReasoning(boolean reasoning) {
		this.reasoning = reasoning;
		}
	
	public boolean isReasoning() {
		return reasoning;
		}	
	
	public boolean isInverseResult() {
		return inverse_result;
		}
	
	public void setInverseResult(boolean inverse_result) {
		this.inverse_result = inverse_result;
		}
	
	@Override
	public void setOutputFile(File outputFile) {
		this.outputFile = outputFile;
		}
	
	public File getOutputFile() {
		return outputFile;
		}
	
	
	@Override
	public int doWork(String[] args)
		{		
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "f:A:Svdo:"))!=-1)
			{
			switch(c)
				{
				case 'o': this.setOutputFile( new File(opt.getOptArg())); break;
				case 'd': this.setReasoning(false); break;
				case 'v': this.inverse_result=true; break;
				case 'S': 
					{
					for(SequenceOntologyTree.Term t:sequenceOntologyTree.getTerms())
						{
						System.out.println(t.getAcn()+"\t"+t.getLabel());
						}
					return 0;
					}
				case 'A':
					{
					this.addTerm(opt.getOptArg().trim());
					break;
					}
				case 'f':
					{
					try
						{
						parseAccessionsFile(new File(opt.getOptArg()));
						}
					catch(IOException err)
						{
						error(err);
						return -1;
						}
					break;
					}
				default: 
					{
					switch(handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		this.initializeKnime();
		List<String> L=new ArrayList<String>();
		for(int i=opt.getOptInd();i<args.length;++i)
			{
			L.add(args[i]);
			}
		return this.executeKnime(L);
		}

	public static void main(String[] args)
		{
		new VcfFilterSequenceOntology().instanceMainWithExit(args);
		}
	}
