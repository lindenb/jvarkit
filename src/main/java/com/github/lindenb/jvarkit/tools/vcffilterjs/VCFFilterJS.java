package com.github.lindenb.jvarkit.tools.vcffilterjs;


/**
 * Author: Pierre Lindenbaum PhD. @yokofakun
 * Motivation http://www.biostars.org/p/66319/ 
 */

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import javax.script.Bindings;
import javax.script.Compilable;
import javax.script.CompiledScript;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;


import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;



public class VCFFilterJS extends AbstractVCFFilter2
	{
	
	
	private CompiledScript  script=null;
	private ScriptEngine engine=null;
	private VCFFilterJS()
		{
		
		}
	
	
	
	
	@Override
	protected void doWork(VcfIterator r, VariantContextWriter w)
			throws IOException
		{
		
		VCFHeader header=r.getHeader();
		
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine( new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		h2.addMetaDataLine( new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
        w.writeHeader(h2);
        Bindings bindings = this.engine.createBindings();
        bindings.put("header", header);
       
        try
	        {
	        while(r.hasNext())
	        	{
	        	VariantContext variation=r.next();
	        	progress.watch(variation.getChr(),variation.getStart());
	        	
				bindings.put("variant", variation);
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
				w.add(variation);
				}
	        }
        catch(ScriptException err)
        	{
        	error(err);
        	throw new IOException(err);
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
				" ( http://sourceforge.net/p/picard/code/HEAD/tree/trunk/src/java/org/broadinstitute/variant/variantcontext/VariantContext.java ) " +
				" and 'header' ( org.broadinstitute.variant.vcf.VCFHeader http://sourceforge.net/p/picard/code/HEAD/tree/trunk/src/java/org/broadinstitute/variant/vcf/VCFHeader.java) in the script context ."
				;
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -e (script) javascript expression.");
		out.println(" -f (script) javascript file.");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		File SCRIPT_FILE=null;
		String SCRIPT_EXPRESSION=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"e:f:"))!=-1)
			{
			switch(c)
				{
				case 'e':SCRIPT_EXPRESSION=opt.getOptArg();break;
				case 'f':SCRIPT_FILE=new File(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		 
		

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
			error("not available: javascript. Do you use the SUN/Oracle JDK ?");
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
			return super.doWork(opt.getOptInd(), args);
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		}
		
	public static void main(String[] args) throws Exception
		{
		new VCFFilterJS().instanceMainWithExit(args);
		}
	
	}
