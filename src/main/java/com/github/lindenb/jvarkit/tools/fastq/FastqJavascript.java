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
* 2014 creation

*/

package com.github.lindenb.jvarkit.tools.fastq;




import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

import javax.script.Bindings;
import javax.script.Compilable;
import javax.script.CompiledScript;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;

import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.CloserUtil;


/**
 * Pierre Lindenbaum PhD @yokofakun
 * @author lindenb
 * filters FASTQ with javascript
 */
public class FastqJavascript
	extends AbstractCommandLineProgram
	{
	private long LIMIT=-1L;

	private CompiledScript  script=null;
	private ScriptEngine engine=null;
	private File failingReadsFile=null;
	private FastqWriter failingReadsWriter=null;
	private boolean interleaved=false;
	
	/** mutable fastq record */
	public static class Record
		{
		private long nLine;
		private String name;
		private String sequence;
		private String name2;
		private String qualities;
		
		Record(FastqRecord rec)
			{
			name=rec.getReadHeader();
			sequence=rec.getReadString();
			name2=rec.getBaseQualityHeader();
			qualities=rec.getBaseQualityString();
			}
		public String getReadHeader() { return name;}
		public String getReadString() { return sequence;}
		public String getBaseQualityHeader() { return name2;}
		public String getBaseQualityString() { return qualities;}
		public long getLine() {
			return nLine;
			}
		public void setName(String name) {
			this.name = name;
			}
		public void setReadString(String sequence) {
			this.sequence = sequence;
		}
		public void setBaseQualityHeader(String name2) {
			this.name2 = name2;
		}
		public void setBaseQualityString(String qualities) {
			this.qualities = qualities;
			}
		private FastqRecord toFastqRecord()
			{
			return new FastqRecord(name,sequence,name2,qualities);
			}
		@Override
		public int hashCode() {
			return toFastqRecord().hashCode();
			}
		@Override
		public String toString() {
			return toFastqRecord().toString();
			}
		}
	
	/** a pair of FASTQ */
	public static class Pair
		{
		private long nLine;
		private Record rec1;
		private Record rec2;
		Pair(Record rec1,Record rec2)
			{
			this.rec1=rec1;
			this.rec2=rec2;
			}
		public long getLine() {
			return nLine;
			}
		public Record getFirst() { return get(0);}
		public Record getSecond() { return get(1);}
		public Record get(int index)
			{
			switch(index)
				{
				case 0: return rec1;
				case 1: return rec2;
				default: throw new IllegalArgumentException("index:="+index);
				}
			}
		
		@Override
		public int hashCode() {
			return get(0).hashCode()*31+get(1).hashCode();
			}
		
		@Override
		public String toString() {
			return get(0).toString()+"\n"+get(1).toString();
			}

		}
	
	private FastqJavascript()
		{
		
		}
	
	/* open failing bam if it was not already open */
	private boolean openFailing()
		{
		if(this.failingReadsFile==null) return false;
		if(this.failingReadsWriter==null)
			{
			info("Writing failings to "+ this.failingReadsFile);
			failingReadsWriter=new BasicFastqWriter(this.failingReadsFile);
			}
		return true;
		}
	
	private void failing(Record rec)
		{
		if(openFailing()) failingReadsWriter.write(rec.toFastqRecord());
		}
	
	private boolean accept(Object result)
		{
		if(result==null)
			{
			return false;
			}
		
		if(result instanceof Boolean)
			{
			if(Boolean.FALSE.equals(result))
				{
				return false;
				}
			}
		else if(result instanceof Number)
			{
			if(((Number)result).intValue()!=1)
				{
				return false;
				}
			}
		else
			{
			warning("script returned neither a number or a boolean "+result.getClass());
			return false;
			}
		return true;
		}
	
	private int doWork(FastqReader r,FastqWriter w)
		{
		try
			{
		
			long count=0L;
	        Bindings bindings = this.engine.createBindings();
	        
			while(r.hasNext())
				{
				Record record=new Record(r.next());
				record.nLine=count;
				
				if(interleaved)
					{
					if(!r.hasNext()) throw new IOException("interleaved: mate missing");
					Record mate= new Record(r.next());
					mate.nLine=count;
					Pair pair=new Pair(record, mate);
					pair.nLine=count;
					bindings.put("pair", pair);
					if(!accept(script.eval(bindings)))
						{
						failing(pair.get(0));
						failing(pair.get(1));
						}
					else
						{
						w.write(pair.get(0).toFastqRecord());
						w.write(pair.get(1).toFastqRecord());
						}
					}
				else
					{
					bindings.put("rec", record);
					if(!accept(script.eval(bindings)))
						{
						failing(record);
						}
					else
						{
						w.write(record.toFastqRecord());
						}
					}
				++count;
				
				if(this.LIMIT>0L && count>=this.LIMIT) break;
				if(System.out.checkError()) break;
				}
			
			openFailing();

			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		}
	
	@Override
	public String getProgramDescription() {
		return "Filters a FASTQ file using javascript( java rhino engine)." +
				"The script puts 'rec' a FastqRecord, or 'pair' for an interleaved input, into the script context .";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/FastqJS";
		}
	
	@Override
	public void printOptions(PrintStream out) {
		out.println(" -f (script) script file");
		out.println(" -e (script) script expression");
		out.println(" -N (limit:int) limit to 'N' records");
		out.println(" -X (fail.fastq) Save dicarded reads in that file. Optional. Default: no file.");
		out.println(" -i interleaved input.");

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
		while((c=opt.getopt(args,getGetOptDefault()+"e:f:N:X:i"))!=-1)
			{
			switch(c)
				{
				case 'X':
					{
					this.failingReadsFile=new File(opt.getOptArg());
					break;
					}
				case 'i': this.interleaved=true;break;
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
		
		FastqWriter fqw=null;
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
			fqw = new BasicFastqWriter(System.out);
			if(opt.getOptInd()==args.length)
				{
				FastqReader in=new FourLinesFastqReader(System.in);
				doWork(in,fqw);
				in.close();
				}
			else for(int i=opt.getOptInd();i< args.length;++i)
				{
				FastqReader in=new FourLinesFastqReader(
						IOUtils.openURIForReading(args[i])
						);
				int ret=doWork(in,fqw);
				in.close();
				if(ret!=0) return ret;
				if(System.out.checkError()) break;
				}
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(failingReadsWriter);
			CloserUtil.close(fqw);
			}
		}
	
	

		
	public static void main(String[] args) throws Exception
		{
		new FastqJavascript().instanceMainWithExit(args);
		}
	
	}
