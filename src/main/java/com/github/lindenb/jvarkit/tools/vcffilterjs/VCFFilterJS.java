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
import java.util.ArrayList;
import java.util.List;

import javax.script.Bindings;
import javax.script.CompiledScript;

import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfTools;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;
import com.google.gson.JsonElement;
import com.google.gson.JsonParser;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;



/**
 * Author: Pierre Lindenbaum PhD. @yokofakun
 * Motivation http://www.biostars.org/p/66319/ 
 */
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

@Program(name="vcffilterjs",description="Filtering VCF with javascript (java Nashorn).")
public class VCFFilterJS
	extends Launcher
	{	
	private static final Logger LOG = Logger.build(VCFFilterJS.class).make();


	
	private CompiledScript compiledScript = null;
	
	
	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;


	@Parameter(names={"-F","--filter"},description="If not empty, variants won't be discarded and this name will be used in the FILTER column")
	private String filteredTag = "";

	@Parameter(names={"-vep","--vep"},description="Decode prediction tag of ensembl VEP")
	private boolean use_vep = false;

	@Parameter(names={"-snpeff","--snpeff"},description="Decode prediction tag of SNPEFF")
	private boolean use_snpeff = false;

	@Parameter(names={"-casecontrol","--casecontrol"},description="Decode Case-Control injected with VcfInjectPedigree")
	private boolean use_casecontrol = false;

	
	@Parameter(names={"-e","--expression"},description=" (js expression). Optional.")
	private String scriptExpr=null;
	@Parameter(names={"-f","--script"},description=" (js file). Optional.")
	private File scriptFile=null;
	@Parameter(names={"-json","--json"},description="json files. syntax key=path/to/file.json . Inject the json object parsed with google gson into the javascript context as 'key'")
	private List<String> jsonFiles=new  ArrayList<>();

	
	public VCFFilterJS()
		{
		
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator r, VariantContextWriter w) {
		try
			{
			
			final VCFHeader header = r.getHeader();
			final VcfTools vcfTools = new VcfTools(header);
			final SnpEffPredictionParser snpEffPredictionParser =(this.use_snpeff ?
					new SnpEffPredictionParserFactory().header(header).get():
					null
					);
			final VepPredictionParser vepPredictionParser = (this.use_vep?
					new VepPredictionParserFactory().header(header).get():
					null
					);

			final  VCFHeader h2 = new VCFHeader(header);
			addMetaData(h2);
			
			final List<Pedigree.Person> individuals =(this.use_casecontrol?
						new ArrayList<>(this.getCasesControlsInPedigree(header)):
						null
						);
			
			final VCFFilterHeaderLine filterHeaderLine = (filteredTag.trim().isEmpty()?null:
				new VCFFilterHeaderLine(this.filteredTag.trim(),"Filtered with "+getProgramName())
				);
			
			
			if(filterHeaderLine!=null) h2.addMetaDataLine(filterHeaderLine);
			
			final  SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header);

			final  Bindings bindings = this.compiledScript.getEngine()
					.createBindings();
			bindings.put("header", header);
			bindings.put("tools", vcfTools);
			
			if(this.use_casecontrol) {
				bindings.put("individuals", individuals);
			}
			
			for(final String jsonkv :this.jsonFiles)
				{
				int eq=jsonkv.indexOf("=");
				if(eq<=0) throw new JvarkitException.UserError("Bad format for json . expected key=/path/to/file.json but got '"+jsonkv+"'");
				final String key=jsonkv.substring(0,eq);
				final FileReader jsonFile = new FileReader(jsonkv.substring(eq+1));
				JsonParser jsonParser=new JsonParser();
				final JsonElement root=jsonParser.parse(jsonFile);
				jsonFile.close();
				bindings.put(key, root);
				}

			w.writeHeader(h2);
			while (r.hasNext() && !w.checkError())
				{
				final  VariantContext variation = progress.watch(r.next());

				bindings.put("variant", variation);
				if(this.use_snpeff)
					{
					bindings.put("snpEff",
						snpEffPredictionParser.getPredictions(variation));
					}
				if(this.use_vep) {
					bindings.put("vep",
						vepPredictionParser.getPredictions(variation));
					}

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
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{

			}
		}
	
	
	@Override
	public int doWork(List<String> args) {
		try 
			{
			this.compiledScript = super.compileJavascript(this.scriptExpr,this.scriptFile);
			
			return doVcfToVcf(args, outputFile);
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			this.compiledScript=null;
			}
		}
	
	
	public static void main(String[] args) throws Exception
		{
		new VCFFilterJS().instanceMainWithExit(args);
		}

	}
