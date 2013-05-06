package fr.inserm.umr1087.jvarkit.tools.samjs;


/**
 * Author: Pierre Lindenbaum PhD. @yokofakun
 * Motivation http://www.biostars.org/p/66319/ 
 */

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;

import javax.script.Bindings;
import javax.script.Compilable;
import javax.script.CompiledScript;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;



public class SamJavascript
	{
	private CompiledScript  script=null;
	private ScriptEngine engine=null;
	private boolean samoutput=false;
	private long limit=-1L;
	private File filenameout=null;
	private SamJavascript()
		{
		
		}
	
	private void scan(InputStream in) throws Exception
		{
		
		long count=0;
		SAMFileReader samFileReader=new SAMFileReader(in);
		samFileReader.setValidationStringency(ValidationStringency.SILENT);
		SAMFileHeader header=samFileReader.getFileHeader();
		
		
        SAMFileWriterFactory sf=new SAMFileWriterFactory();
        sf.setCreateIndex(false);
        File stdout=(filenameout==null?new File("/dev/stdout"):filenameout);
     	
        SAMFileWriter sw;
        if(!samoutput )
        	{
        	if(!stdout.exists() && filenameout==null) //stdout
			{
			throw new IOException("Cannot save as BAM because "+stdout+" doesn't exist. Please use SAM.");
			}
        	sw=sf.makeBAMWriter(header,false,stdout);        
        	}
        else if(filenameout==null)
        	{
        	sw=sf.makeSAMWriter(header,false,System.out);        
        	}
	else
		{
		sw=sf.makeSAMWriter(header,false,filenameout);
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
			else if(result instanceof Integer)
				{
				if(((Integer)result).intValue()!=1) continue;
				}
			else
				{
				continue;
				}
			++count;
			sw.addAlignment(record);
			if(this.limit!=-1L && count>=this.limit) break;
			}
		samFileReader.close();
		sw.close();
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
				System.err.println(" -S : SAM output.");
				System.err.println(" -e (script).");
				System.err.println(" -o (fileout.bam) or stdout");
				System.err.println(" -f (file).");
				System.err.println(" -L (int) limit to 'L' records.");
				System.err.println("the script puts 'record' a SamRecord (http://picard.sourceforge.net/javadoc/net/sf/samtools/SAMRecord.html)  " +
						" and 'header' ( http://picard.sourceforge.net/javadoc/net/sf/samtools/SAMFileHeader.html) in the script context .");
				return;
				}
			else if(args[optind].equals("-S") )
				{
				this.samoutput=true;
				}
			else if(args[optind].equals("-L") && optind+1< args.length)
				{
				this.limit=Long.parseLong(args[++optind]);
				}
			else if(args[optind].equals("-o") && optind+1< args.length)
				{
				this.filenameout=new File(args[++optind]);
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
			scan(System.in);
			}
		else if(optind+1==args.length)
			{
			FileInputStream fin=new FileInputStream(args[optind]);
			scan(fin);
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
		new SamJavascript().run(args);
		}
	
	}
