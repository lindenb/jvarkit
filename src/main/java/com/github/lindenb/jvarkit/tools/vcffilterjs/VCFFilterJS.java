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

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;

import org.broad.tribble.readers.LineReader;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.tools.vcffixindels.VCFFixIndels;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter;


public class VCFFilterJS extends AbstractVCFFilter
	{
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+
	" Filtering VCF with javascript (java rhino)."+
	" The script puts 'variant' a org.broadinstitute.variant.variantcontext.VariantContext " +
	" ( http://sourceforge.net/p/picard/code/HEAD/tree/trunk/src/java/org/broadinstitute/variant/variantcontext/VariantContext.java ) " +
	" and 'header' ( org.broadinstitute.variant.vcf.VCFHeader http://sourceforge.net/p/picard/code/HEAD/tree/trunk/src/java/org/broadinstitute/variant/vcf/VCFHeader.java) in the script context ."
	;
	@Option(shortName="SF", doc="javascript file ",optional=true)
	public File SCRIPT_FILE=null;
	@Option(shortName="SE", doc="javascript expression ",optional=true)
	public String SCRIPT_EXPRESSION=null;
	
	private static final Log LOG=Log.getInstance(VCFFixIndels.class);

	
	private CompiledScript  script=null;
	private ScriptEngine engine=null;
	public VCFFilterJS()
		{
		
		}
	
	@Override
	public String getVersion() {
		return "1.0";
		}
	
	
	@Override
	protected void doWork(LineReader in, VariantContextWriter w)
			throws IOException
		{
		
		VCFCodec vcfCodec=new VCFCodec();
		VCFHeader header=(VCFHeader)vcfCodec.readHeader(in);
		
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getName(),"Filtered with "+(SCRIPT_FILE==null?"":SCRIPT_FILE)+(SCRIPT_EXPRESSION==null?"":SCRIPT_EXPRESSION)));
		
		
        w.writeHeader(h2);
        Bindings bindings = this.engine.createBindings();
        bindings.put("header", header);
       
        try
	        {
	        String line;
	        while((line=in.readLine())!=null)
	        	{
	        	VariantContext variation=vcfCodec.decode(line);
	        	
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
					continue;
					}
				w.add(variation);
				}
	        }
        catch(ScriptException err)
        	{
        	LOG.error(err);
        	throw new IOException(err);
        	}
		}

	@Override
	protected int doWork()
		{
		try
			{
			if(SCRIPT_EXPRESSION==null && SCRIPT_FILE==null)
				{
				LOG.error("undefined script");
				return -1;
				}
			if(SCRIPT_EXPRESSION!=null && SCRIPT_FILE!=null)
				{
				LOG.error("both file/expr are set");
				return -1;
				}
			ScriptEngineManager manager = new ScriptEngineManager();
			this.engine = manager.getEngineByName("js");
			if(this.engine==null)
				{
				LOG.error("not available: javascript. Use the SUN/Oracle JDK ?");
				return -1;
				}
			
			Compilable compilingEngine = (Compilable)this.engine;
			this.script = null;
			if(SCRIPT_FILE!=null)
				{
				LOG.info("Compiling "+SCRIPT_FILE);
				FileReader r=new FileReader(SCRIPT_FILE);
				this.script=compilingEngine.compile(r);
				r.close();
				}
			else
				{
				LOG.info("Compiling "+SCRIPT_EXPRESSION);
				this.script=compilingEngine.compile(SCRIPT_EXPRESSION);
				}
			}
		catch(Exception err)
			{
			LOG.error(err,"Error occured");
			}
		return super.doWork();
		}
		
	public static void main(String[] args) throws Exception
		{
		new VCFFilterJS().instanceMainWithExit(args);
		}
	
	}
