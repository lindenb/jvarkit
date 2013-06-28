package com.github.lindenb.jvarkit.tools.samjs;


/**
 * Author: Pierre Lindenbaum PhD. @yokofakun
 * Motivation http://www.biostars.org/p/66319/ 
 */

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;

import javax.script.Bindings;
import javax.script.Compilable;
import javax.script.CompiledScript;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;


import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;



public class SamJavascript extends CommandLineProgram
	{
	private static final Log LOG=Log.getInstance(SamJavascript.class);
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"Filters a BAM using javascript( java rhino engine)." +
			"The script puts 'record' a SamRecord (http://picard.sourceforge.net/javadoc/net/sf/samtools/SAMRecord.html)  " +
			" and 'header' ( http://picard.sourceforge.net/javadoc/net/sf/samtools/SAMFileHeader.html) in the script context .";
    
	@Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="BAM file to process. Default stdin. ",optional=true)
	public File IN=null;
	@Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="output filename. Default stdout. ",optional=true)
	public File OUT=null;
	@Option(shortName="SF", doc="javascript file ",optional=true)
	public File SCRIPT_FILE=null;
	@Option(shortName="SE", doc="javascript expression ",optional=true)
	public String SCRIPT_EXPRESSION=null;
	@Option(shortName="SAM",doc="sam output ",optional=true)
	public boolean SAM_OUTPUT=false;
	@Option(shortName="L",doc="limit to 'L' records. ",optional=true)
	public Long LIMIT=null;
	
	
	
	private CompiledScript  script=null;
	private ScriptEngine engine=null;
	
	@Override
	public String getVersion() {
		return "1.0";
		}
	
	private void scan(SAMFileReader samFileReader) throws Exception
		{
		
		long count=0;
		
		samFileReader.setValidationStringency(super.VALIDATION_STRINGENCY);
		SAMFileHeader header=samFileReader.getFileHeader();
		
		
        SAMFileWriterFactory sf=new SAMFileWriterFactory();
        sf.setCreateIndex(false);
        File stdout=(OUT==null?new File("/dev/stdout"):OUT);
     	
        SAMFileWriter sw=null;
        if(!SAM_OUTPUT )
        	{
        	if(!stdout.exists() && OUT==null) //stdout
				{
        		samFileReader.close();
				throw new IOException("Cannot save as BAM because "+stdout+" doesn't exist. Please use SAM.");
				}
        	sw=sf.makeBAMWriter(header,false,stdout);        
        	}
        else if(OUT==null)
        	{
        	sw=sf.makeSAMWriter(header,false,System.out);        
        	}
		else
			{
			sw=sf.makeSAMWriter(header,false,OUT);
			}

        Bindings bindings = this.engine.createBindings();
        bindings.put("header", header);
       
        
		for(Iterator<SAMRecord> iter=samFileReader.iterator();
				iter.hasNext(); )
			{
			SAMRecord record=iter.next();
			bindings.put("record", record);
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
				LOG.info("script returned neither a number or a boolean "+result.getClass());
				continue;
				}
			++count;
			sw.addAlignment(record);
			if(this.LIMIT!=null && count>=this.LIMIT) break;
			}
		sw.close();
		}
	
	@Override
	public int doWork()
		{
		try {
			if(SCRIPT_EXPRESSION==null && SCRIPT_FILE==null)
				{
				LOG.error("undefined script");
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
				FileReader r=new FileReader(SCRIPT_FILE);
				this.script=compilingEngine.compile(r);
				r.close();
				}
			else
				{
				this.script=compilingEngine.compile(SCRIPT_EXPRESSION);
				}
			
			if(IN==null)
				{
				SAMFileReader samFileReader=new SAMFileReader(System.in);
				scan(samFileReader);
				samFileReader.close();
				}
			else
				{
				SAMFileReader samFileReader=new SAMFileReader(IN);
				scan(samFileReader);
				samFileReader.close();
				}
			
			return 0;
			
		} catch (Exception e) {
			LOG.error(e);
			return -1;
			}
		}
	

		
	public static void main(String[] args) throws Exception
		{
		new SamJavascript().instanceMainWithExit(args);
		}
	
	}
