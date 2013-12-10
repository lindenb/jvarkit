package com.github.lindenb.jvarkit.tools.samjs;


/**
 * Author: Pierre Lindenbaum PhD. @yokofakun
 * Motivation http://www.biostars.org/p/66319/ 
 */

import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.Iterator;

import javax.script.Bindings;
import javax.script.Compilable;
import javax.script.CompiledScript;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;
import com.github.lindenb.jvarkit.util.picard.AbstractBamFilterProgram;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;



public class SamJavascript
	extends AbstractBamFilterProgram
	{
	private long LIMIT=-1L;

	private CompiledScript  script=null;
	private ScriptEngine engine=null;
	
	private SamJavascript()
		{
		
		}
	
	@Override
	protected int doWork(SAMFileReader samFileReader, SAMFileWriter sw)
		{
		try
			{
			long count=0L;
	        Bindings bindings = this.engine.createBindings();
	        bindings.put("header", samFileReader.getFileHeader());
	       
	        
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
					warning("script returned neither a number or a boolean "+result.getClass());
					continue;
					}
				++count;
				sw.addAlignment(record);
				if(this.LIMIT>0L && count>=this.LIMIT) break;
				}
			sw.close();
			return 0;
			}
		catch(ScriptException err)
			{
			error(err);
			return -1;
			}
		}
	
	@Override
	public String getProgramDescription() {
		return "Filters a BAM using javascript( java rhino engine)." +
				"The script puts 'record' a SamRecord (http://picard.sourceforge.net/javadoc/net/sf/samtools/SAMRecord.html)  " +
				" and 'header' ( http://picard.sourceforge.net/javadoc/net/sf/samtools/SAMFileHeader.html) in the script context .";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/SamJS";
		}
	
	@Override
	public void printOptions(PrintStream out) {
		out.println(" -f (script) script file");
		out.println(" -e (script) script expression");
		out.println(" -N (limit:int) limit to 'N' records");
		super.printOptions(out);
		}

	
	@Override
	public int doWork(String[] args)
		{
		try
			{
			ScriptEngineManager manager = new ScriptEngineManager();
			this.engine = manager.getEngineByName("js");
			if(this.engine==null)
				{
				error("not available: javascript. Use the SUN/Oracle JDK ?");
				return -1;
				}
		
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
			
		
		String SCRIPT_EXPRESSION=null;
		File SCRIPT_FILE=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"e:f:N:"))!=-1)
			{
			switch(c)
				{
				case 'e':SCRIPT_EXPRESSION=opt.getOptArg();break;
				case 'f':SCRIPT_FILE=new File(opt.getOptArg());break;
				case 'N': this.LIMIT=Long.parseLong(opt.getOptArg());break;
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
		try
			{
			Compilable compilingEngine = (Compilable)this.engine;
			this.script = null;
			if(SCRIPT_FILE!=null)
				{
				info("Reading script "+SCRIPT_FILE);
				FileReader r=new FileReader(SCRIPT_FILE);
				this.script=compilingEngine.compile(r);
				r.close();
				}
			else
				{
				info("Eval script "+SCRIPT_EXPRESSION);
				this.script=compilingEngine.compile(SCRIPT_EXPRESSION);
				}
			return super.doWork(this.bamFileOut, opt.getOptInd(), args);
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		}
	
	

		
	public static void main(String[] args) throws Exception
		{
		new SamJavascript().instanceMainWithExit(args);
		}
	
	}
