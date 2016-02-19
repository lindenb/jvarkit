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


import java.io.IOException;
import java.util.Collection;

import javax.script.Bindings;
import javax.script.CompiledScript;
import javax.script.ScriptException;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;



/**
 * Author: Pierre Lindenbaum PhD. @yokofakun
 * Motivation http://www.biostars.org/p/66319/ 
 */
public class VCFFilterJS
	extends AbstractVCFFilterJS
	{	
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(AbstractVCFFilterJS.class);
	private CompiledScript compiledScript = null;
	
	/** 2015-02-10 : moved to public , so we can use it in knime */
	public VCFFilterJS()
		{
		
		}
	@Override
	public Collection<Throwable> doVcfToVcf(
			final String inputName,
			final VcfIterator r,
			final VariantContextWriter w
			) throws IOException
		{
		try
			{
			final VCFHeader header = r.getHeader();
			final SnpEffPredictionParser snpEffPredictionParser = new SnpEffPredictionParser(
					header);
			final VepPredictionParser vepPredictionParser = new VepPredictionParser(
					header);

			final  VCFHeader h2 = new VCFHeader(header);
			addMetaData(h2);
			
			final VCFFilterHeaderLine filterHeaderLine = (filteredTag.trim().isEmpty()?null:
				new VCFFilterHeaderLine(super.filteredTag.trim(),"Filtered with "+getName())
				);
			
			
			if(filterHeaderLine!=null) h2.addMetaDataLine(filterHeaderLine);
			
			final  SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header);

			final  Bindings bindings = this.compiledScript.getEngine()
					.createBindings();
			bindings.put("header", header);

			w.writeHeader(h2);
			while (r.hasNext() && !w.checkError())
				{
				final  VariantContext variation = progress.watch(r.next());

				bindings.put("variant", variation);
				bindings.put("snpEff",
						snpEffPredictionParser.getPredictions(variation));
				bindings.put("vep",
						vepPredictionParser.getPredictions(variation));

				if (!evalJavaScriptBoolean(this.compiledScript, bindings) )
					{
					if(filterHeaderLine!=null)
						{
						final VariantContextBuilder vcb = new VariantContextBuilder(variation);
						vcb.filter(filterHeaderLine.getID());
						w.add(vcb.make());
						}
					continue;
					}
				w.add(variation);
				}
			return RETURN_OK;
			}
		catch (ScriptException err)
			{
			return wrapException(err);
			}
		finally
			{

			}
		}
	
	@Override
	public Collection<Throwable> initializeKnime()
		{
		try
			{
			this.compiledScript = super.compileJavascript();
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		return super.initializeKnime();
		}
	
	@Override
	public void disposeKnime()
		{
		this.compiledScript=null;
		super.disposeKnime();
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		LOG.info("NPUT ="+inputName+" "+getInputFiles());
		return doVcfToVcf(inputName);
		}
	
	public static void main(String[] args) throws Exception
		{
		new VCFFilterJS().instanceMainWithExit(args);
		}

	}
