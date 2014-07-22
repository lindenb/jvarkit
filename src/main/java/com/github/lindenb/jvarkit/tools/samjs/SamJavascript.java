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
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;



public class SamJavascript
	extends AbstractBamFilterProgram
	{
	private long LIMIT=-1L;

	private CompiledScript  script=null;
	private ScriptEngine engine=null;
	private File failingReadsFile=null;
	private SAMFileWriter failingReadsWriter=null;
	

	private SamJavascript()
		{
		
		}
	
	/* open failing bam if it was not already open */
	private void openFailing(SAMFileHeader h)
		{
		if(this.failingReadsFile==null) return;
		if(this.failingReadsWriter==null)
			{
			info("Writing failings to "+ this.failingReadsFile);
			SAMFileHeader header=createOuputSamFileHeader(h);
			
			SAMFileWriterFactory sfwf=new SAMFileWriterFactory();
			sfwf.setCreateIndex(super.create_index);
			sfwf.setCreateMd5File(super.create_md5);
			sfwf.setMaxRecordsInRam(super.max_records_in_ram);
			

			
			failingReadsWriter=sfwf.makeSAMOrBAMWriter(
					header, 
					isOutputPresorted(header)
					,this.failingReadsFile);
			}
		}
	
	private void failing(SAMRecord rec,SAMFileHeader h)
		{
		openFailing(h);
		if(failingReadsWriter!=null) failingReadsWriter.addAlignment(rec);
		}
	
	@Override
	protected int doWork(SamReader samFileReader, SAMFileWriter sw)
		{
		try
			{
			SAMFileHeader header=samFileReader.getFileHeader();
			long count=0L;
	        Bindings bindings = this.engine.createBindings();
	        bindings.put("header", samFileReader.getFileHeader());
	        SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
	        
			for(Iterator<SAMRecord> iter=samFileReader.iterator();
					iter.hasNext(); )
				{
				SAMRecord record=iter.next();
				progress.watch(record);
				bindings.put("record", record);
				Object result = script.eval(bindings);
				if(result==null)
					{
					failing(record,header);
					continue;
					}
				
				if(result instanceof Boolean)
					{
					if(Boolean.FALSE.equals(result))
						{
						failing(record,header);
						continue;
						}
					}
				else if(result instanceof Number)
					{
					if(((Number)result).intValue()!=1)
						{
						failing(record,header);
						continue;
						}
					}
				else
					{
					warning("script returned neither a number or a boolean "+result.getClass());
					failing(record,header);
					continue;
					}
				++count;
				sw.addAlignment(record);
				if(this.LIMIT>0L && count>=this.LIMIT) break;
				}
			sw.close();
			
			openFailing(header);
			
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
				"The script puts 'record' a SamRecord (http://picard.sourceforge.net/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html)  " +
				" and 'header' ( http://picard.sourceforge.net/javadoc/htsjdk/htsjdk/samtools/SAMFileHeader.html ) in the script context .";
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
		out.println(" -X (fail.bam) Save dicarded reads in that file. Optional. Default: no file.");

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
		while((c=opt.getopt(args,getGetOptDefault()+"e:f:N:X:"))!=-1)
			{
			switch(c)
				{
				case 'X':
					{
					this.failingReadsFile=new File(opt.getOptArg());
					break;
					}
				case 'e':SCRIPT_EXPRESSION=opt.getOptArg();break;
				case 'f':SCRIPT_FILE=new File(opt.getOptArg());break;
				case 'N': this.LIMIT=Long.parseLong(opt.getOptArg());break;
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
		finally
			{
			CloserUtil.close(failingReadsWriter);
			}
		}
	
	

		
	public static void main(String[] args) throws Exception
		{
		new SamJavascript().instanceMainWithExit(args);
		}
	
	}
