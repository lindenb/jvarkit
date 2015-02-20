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

* 2015 moving to knime and adding predictions snpeff and VEP
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.vcffilterjs;



import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.script.Bindings;
import javax.script.Compilable;
import javax.script.CompiledScript;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;

import com.github.lindenb.jvarkit.knime.KnimeApplication;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;



/**
 * Author: Pierre Lindenbaum PhD. @yokofakun
 * Motivation http://www.biostars.org/p/66319/ 
 */
public class VCFFilterJS
	extends AbstractCommandLineProgram
	implements KnimeApplication
	{
	private CompiledScript  script=null;
	private ScriptEngine engine=null;
	/** output file : stdout if null */
	private File outputFile=null;
	/** number of variants filterer */
	private int countFilteredVariants=0;
	/** expression in file */
	private File SCRIPT_FILE=null;
	/** expression in string */
	private String SCRIPT_EXPRESSION=null;

	
	/** 2015-02-10 : moved to public , so we can use it in knime */
	public VCFFilterJS()
		{
		
		}
	
	public void setScriptExpression(String expression) {
		SCRIPT_EXPRESSION = expression;
		}
	
	public void setScriptFile(File src)
		{
		SCRIPT_FILE = src;
		}
	
	/** for knime, return the number of variants kept after execute*/ 
	public int getVariantCount()
		{
		return this.countFilteredVariants;
		}

	private void filterVcfIterator(VcfIterator r) throws IOException
		{
		this.countFilteredVariants=0;
		VariantContextWriter w = null;
		try
			{
			VCFHeader header=r.getHeader();
			SnpEffPredictionParser snpEffPredictionParser=new SnpEffPredictionParser(header);
			VepPredictionParser vepPredictionParser=new VepPredictionParser(header);
			
			
			VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
	
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
	
			
			if(getOutputFile()==null)
				{
				w = VCFUtils.createVariantContextWriterToStdout();
				}
			else
				{
				info("opening vcf writer to "+getOutputFile());
				w = VCFUtils.createVariantContextWriter(getOutputFile());
				}
	
			
			Bindings bindings = this.engine.createBindings();
	        bindings.put("header", header);
       
	        w.writeHeader(h2);
	        while(r.hasNext())
	        	{
	        	VariantContext variation=progress.watch(r.next());
	        	
				bindings.put("variant", variation);
				
				
				bindings.put("snpEff", snpEffPredictionParser.getPredictions(variation));
				bindings.put("vep", vepPredictionParser.getPredictions(variation));
				
				Object result = script.eval(bindings);
				if(result==null) continue;
				if(result instanceof Boolean)
					{
					if(Boolean.FALSE.equals(result)) continue;
					}
				else if(result instanceof Number)
					{
					if(((Number)result).intValue()!=1) continue;
					}
				else
					{
					warning("Script returned something that is not a boolean or a number:"+result.getClass());
					continue;
					}
				this.countFilteredVariants++;
				w.add(variation);
				}
	        }
        catch(ScriptException err)
        	{
        	error(err);
        	throw new IOException(err);
        	}
		finally
			{
			CloserUtil.close(w);
			}
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VCFFilterJS";
		}
	
	@Override
	public String getProgramDescription() {
		return  "Filtering VCF with javascript (java rhino)."+
				" The script puts 'variant' a org.broadinstitute.variant.variantcontext.VariantContext " +
				" ( https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html ) " +
				" , 'header' ( org.broadinstitute.variant.vcf.VCFHeader https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/vcf/VCFHeader.html ) in the script context, "+
				" 'snpEff' a java.util.List SnpEffPredictionParser$SnpEffPrediction of https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/util/vcf/predictions/SnpEffPredictionParser.java "+
				" 'vep' a java.util.List of VepPredictionParser$VepPrediction https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/util/vcf/predictions/VepPredictionParser.java  ."
				;
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -e (script) javascript expression.");
		out.println(" -f (script) javascript file.");
		out.println(" -o (file) output file. default: stdout");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"e:f:o:"))!=-1)
			{
			switch(c)
				{
				case 'o':this.setOutputFile(new File(opt.getOptArg()));break;
				case 'e':this.setScriptExpression(opt.getOptArg());break;
				case 'f':this.setScriptFile(new File(opt.getOptArg()));break;
				default:
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		try
			{
			if(this.initializeKnime()!=0)
				{
				return -1;
				}
			List<String> L=new ArrayList<String>();
			for(int i=opt.getOptInd();i<args.length;++i)
				{
				L.add(args[i]);
				}
			return this.executeKnime(L);
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		}
		


	@Override
	public int initializeKnime()
		{
		if(SCRIPT_EXPRESSION==null && SCRIPT_FILE==null)
			{
			error("undefined script");
			return -1;
			}
		if(SCRIPT_EXPRESSION!=null && SCRIPT_FILE!=null)
			{
			error("both file/expr are set");
			return -1;
			}
		ScriptEngineManager manager = new ScriptEngineManager();
		this.engine = manager.getEngineByName("js");
		if(this.engine==null)
			{
			error("The embedded 'javascript' engine is not available in java. Do you use the SUN/Oracle Java Runtime ?");
			return -1;
			}
		
		try
			{
			Compilable compilingEngine = (Compilable)this.engine;
			this.script = null;
			if(SCRIPT_FILE!=null)
				{
				info("Compiling "+SCRIPT_FILE);
				FileReader r=new FileReader(SCRIPT_FILE);
				this.script=compilingEngine.compile(r);
				r.close();
				}
			else
				{
				info("Compiling "+SCRIPT_EXPRESSION);
				this.script=compilingEngine.compile(SCRIPT_EXPRESSION);
				}
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		return 0;
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


	@Override
	public void disposeKnime() {
		this.engine=null;
		this.SCRIPT_EXPRESSION=null;
		this.SCRIPT_FILE=null;
		
	}


	@Override
	public void checkKnimeCancelled() {
		}


	@Override
	public void setOutputFile(File out) {
		this.outputFile=out;
		}

	public File getOutputFile() {
		return outputFile;
		}
	
	public static void main(String[] args) throws Exception
		{
		new VCFFilterJS().instanceMainWithExit(args);
		}

	}
