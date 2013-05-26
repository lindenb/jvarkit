package com.github.lindenb.jvarkit.tools.vcffilterjs;


/**
 * Author: Pierre Lindenbaum PhD. @yokofakun
 * Motivation http://www.biostars.org/p/66319/ 
 */

import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.util.EnumSet;

import javax.script.Bindings;
import javax.script.Compilable;
import javax.script.CompiledScript;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;

import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextUtils;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;


public class VCFFilterJS
	{
	private CompiledScript  script=null;
	private ScriptEngine engine=null;
	private VCFFilterJS()
		{
		
		}
	
	private void scan(LineReader in) throws Exception
		{
		
		VCFCodec vcfCodec=new VCFCodec();
		VCFHeader header=(VCFHeader)vcfCodec.readHeader(in);
		EnumSet<Options> options=EnumSet.noneOf(Options.class);
		VariantContextWriter w=VariantContextWriterFactory.create(
				System.out,
				null,
				options
				);
		
        w.writeHeader(header);
        Bindings bindings = this.engine.createBindings();
        bindings.put("header", header);
       
        
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
	
	
	private void run(String[] args)
		throws Exception
		{
		String scriptStr=null;
		File scriptFile=null;
		int optind=0;
		while(optind< args.length)
			{
			if(args[optind].equals("-h") ||
			   args[optind].equals("-help") ||
			   args[optind].equals("--help"))
				{
				System.err.println("Pierre Lindenbaum PhD. 2013");
				System.err.println("Options:");
				System.err.println(" -h help; This screen.");
				System.err.println(" -e (script).");
				System.err.println(" -f (file).");
				System.err.println(
					"the script puts 'variant' a org.broadinstitute.variant.variantcontext.VariantContext " +
						" ( http://sourceforge.net/p/picard/code/HEAD/tree/trunk/src/java/org/broadinstitute/variant/variantcontext/VariantContext.java ) " +
						" and 'header' ( org.broadinstitute.variant.vcf.VCFHeader http://sourceforge.net/p/picard/code/HEAD/tree/trunk/src/java/org/broadinstitute/variant/vcf/VCFHeader.java) in the script context .");
				return;
				}
			else if(args[optind].equals("-e") && optind+1< args.length)
				{
				scriptStr=args[++optind];
				}
			else if(args[optind].equals("-f") && optind+1< args.length)
				{
				scriptFile=new File(args[++optind]);
				}
			else if(args[optind].equals("--"))
				{
				optind++;
				break;
				}
			else if(args[optind].startsWith("-"))
				{
				System.err.println("Unknown option "+args[optind]);
				return;
				}
			else 
				{
				break;
				}
			++optind;
			}
		if(scriptStr==null && scriptFile==null)
			{
			System.err.println("undefined script");
			System.exit(-1);
			}
		
		ScriptEngineManager manager = new ScriptEngineManager();
		this.engine = manager.getEngineByName("js");
		if(this.engine==null)
			{
			System.err.println("not available: javascript. Use the SUN/Oracle JDK ?");
			System.exit(-1);
			}
		
		Compilable compilingEngine = (Compilable)this.engine;
		this.script = null;
		if(scriptFile!=null)
			{
			FileReader r=new FileReader(scriptFile);
			this.script=compilingEngine.compile(r);
			r.close();
			}
		else
			{
			this.script=compilingEngine.compile(scriptStr);
			}
		
		if(optind==args.length)
			{
			scan(new AsciiLineReader(System.in));
			}
		else if(optind+1==args.length)
			{
			FileInputStream fin=new FileInputStream(args[optind]);
			scan(new AsciiLineReader(fin));
			fin.close();
			}
		else 
			{
			System.err.println("illegal number of arguments.");
			System.exit(-1);
			}
		}
		
	public static void main(String[] args) throws Exception
		{
		new VCFFilterJS().run(args);
		}
	
	}
