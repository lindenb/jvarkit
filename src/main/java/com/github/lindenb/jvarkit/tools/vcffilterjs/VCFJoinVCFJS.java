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
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

import javax.script.Bindings;
import javax.script.ScriptException;

import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.IndexedVcfFileReader;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;



/**
 * Author: Pierre Lindenbaum PhD. @yokofakun
 */
public class VCFJoinVCFJS
	extends AbstractVcfJavascript
	{	
	private File databaseFile=null;
	private IndexedVcfFileReader indexedVcfFileReader=null;

	public VCFJoinVCFJS()
		{
		
		}
	
	public void setDatabaseFile(File databaseFile) {
		this.databaseFile = databaseFile;
		}
	
	
	@Override
	public int initializeKnime() {
		if(databaseFile==null)
			{
			error("Database file undefined");
			return -1;
			}
		if(!databaseFile.exists() || !databaseFile.isFile() || !databaseFile.canRead())
			{
			error("Cannot read VCF file "+this.databaseFile);
			return -1;
			}
		try
			{
			info("Opening datbase file "+this.databaseFile);
			indexedVcfFileReader =new IndexedVcfFileReader(this.databaseFile);
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		return super.initializeKnime();
		}
	
	@Override
	public void disposeKnime() {
		CloserUtil.close(this.indexedVcfFileReader);
		this.indexedVcfFileReader=null;
		super.disposeKnime();
		}
	
	@Override
	protected void doWork(
			String inpuSource,
			VcfIterator r,
			VariantContextWriter w
			) throws IOException
		{
		try
			{
			
			VCFHeader header=r.getHeader();
			
			
			VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
			
			
			SAMSequenceDictionary dict1=header.getSequenceDictionary();
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict1);
				
			Bindings bindings = this.getScriptEngine().createBindings();
	        bindings.put("header", header);
	        
	        
	        bindings.put("dbheader", indexedVcfFileReader.getHeader());
			SAMSequenceDictionary dict2=indexedVcfFileReader.getHeader().getSequenceDictionary();
			if(dict1!=null && dict2!=null)
				{
				boolean found=false;
				for(SAMSequenceRecord ssr:dict1.getSequences())
					{
					if(dict2.getSequence(ssr.getSequenceName())==null) continue;
					found=true;
					break;
					}
				if(!found)
					{
					warning("WARNING : no common chromosomes between "+databaseFile+" and "+inpuSource);
					}
				}
			
			
	        
	        w.writeHeader(h2);
	        while(r.hasNext() && !this.checkOutputError())
	        	{
	        	VariantContext variation=progress.watch(r.next());
	        	
	        	List<VariantContext> joinVariant = indexedVcfFileReader.getVariants(
	        			variation.getContig(),
	        			variation.getStart(),
	        			variation.getEnd()
	        			);
				bindings.put("variant", variation);
				bindings.put("dbvariants", joinVariant);
				
				if(!accept(bindings)) continue;				

				this.incrVariantCount();
				w.add(variation);
				
				if(super.checkOutputError()) break;
				}
	        }
        catch(ScriptException err)
        	{
        	error(err);
        	throw new IOException(err);
        	}
		finally
			{
			}
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"VCFJoinVcfJS";
		}
	
	@Override
	public String getProgramDescription() {
		return  "Join two VCF files and accept/exclude variant with javascript (java rhino)."+
				" The script injects in the javascript context :\n"+
				"* 'variant' "+VariantContext.class.getName()+" ("+HtsjdkVersion.getJavadocUrl(VariantContext.class)+")\n"+
				"* 'header' ( "+VCFHeader.class.getName()+" "+ HtsjdkVersion.getJavadocUrl(VCFHeader.class)+" )\n"+
				"* 'dbvariants' java.util.List<"+VariantContext.class.getName()+"> ("+HtsjdkVersion.getJavadocUrl(VariantContext.class)+")\n"+
				"* 'dbheader' ( "+VCFHeader.class.getName()+" "+ HtsjdkVersion.getJavadocUrl(VCFHeader.class)+" )"
				;
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -d (vcf database file) a tribble or a tabix file. REQUIRED." );
		super.printOptions(out);
		}
	
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"d:"))!=-1)
			{
			switch(c)
				{
				case 'd': setDatabaseFile(new File(opt.getOptArg())); break;
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
		return mainWork(opt.getOptInd(), args);
		}
		

	
	public static void main(String[] args) throws Exception
		{
		new VCFJoinVCFJS().instanceMainWithExit(args);
		}

	}
