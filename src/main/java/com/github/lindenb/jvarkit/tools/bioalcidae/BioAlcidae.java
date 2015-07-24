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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.bioalcidae;


import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.readers.LineIterator;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;

import javax.script.Bindings;
import javax.script.Compilable;
import javax.script.CompiledScript;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.knime.AbstractKnimeApplication;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;



/**
 * BioAlcidae
 *
 */
public class BioAlcidae
	extends AbstractKnimeApplication
	{
	private enum FORMAT {
		VCF{
			@Override
			boolean canAs(String src) {
				return src!=null && (src.endsWith(".vcf") || src.endsWith(".vcf.gz") );
			}
			},
		SAM{
			@Override
			boolean canAs(String src) {
				return src!=null && (src.endsWith(".sam") || src.endsWith(".bam") );
			}
			},
		BAM{
			@Override
			boolean canAs(String src) {
				return src!=null && (src.endsWith(".sam") || src.endsWith(".bam") );
			}
			},
		FASTA{
			@Override
			boolean canAs(String src) {
				return src!=null && (src.endsWith(".fa") || src.endsWith(".fasta") || src.endsWith(".fa.gz") || src.endsWith(".fasta.gz")  );
			}
			},
		FASTQ{
				@Override
				boolean canAs(String src) {
					return src!=null && (src.endsWith(".fq") || src.endsWith(".fastq") || src.endsWith(".fq.gz") || src.endsWith(".fastq.gz")  );
				}
				}
			
			;
		abstract boolean canAs(String src);
		};
	private Bindings bindings=null;
	private CompiledScript  script=null;
	private ScriptEngine engine=null;
	private FORMAT format= null;
	private PrintWriter writer=null;

	/** expression in file */
	private File SCRIPT_FILE=null;
	/** expression in string */
	private String SCRIPT_EXPRESSION=null;	
	
	
	public BioAlcidae()
		{
		
		}
	@Override
	public String getProgramDescription() {
		return "javascript version of awk";
		}			
	
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"BioAlcidae";
		}
	
	
	public void setScriptExpression(String expression) {
		SCRIPT_EXPRESSION = expression;
		}
	
	public void setScriptFile(File src)
		{
		SCRIPT_FILE = src;
		}
	
	public void setFormat(String fmt )
		{
		try {
			this.format=FORMAT.valueOf(fmt.toUpperCase());
		} catch (Exception e) {
			this.format=null;
			}
		}

	@Override
	public int initializeKnime()
		{
		if(SCRIPT_EXPRESSION==null && SCRIPT_FILE==null)
			{
			error("undefined script");
			return -1;
			}
		if(SCRIPT_EXPRESSION!=null && SCRIPT_FILE!=null)
			{
			error("both javascript file/expr are set");
			return -1;
			}
		ScriptEngineManager manager = new ScriptEngineManager();
		this.engine = manager.getEngineByName("js");
		if(this.engine==null)
			{
			error("The embedded 'javascript' engine is not available in java. Do you use the SUN/Oracle Java Runtime ?");
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
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		return 0;
		}
	
	private int execute_vcf(String source) throws IOException
		{
		info("source: "+source);
		VcfIterator in=null;
		try {
			in = VCFUtils.createVcfIterator(source);
			bindings.put("codec",in.getCodec());
			bindings.put("header",in.getHeader());
			bindings.put("iter",in);
			bindings.put("format","vcf");
			this.script.eval(bindings);
			return 0;
			} 
		catch (Exception e)
			{
			error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			bindings.remove("format");
			bindings.remove("codec");
			bindings.remove("header");
			bindings.remove("iter");
			}
		}
	
	private int execute_bam(String source) throws IOException
		{
		SamReader in=null;
		SAMRecordIterator iter=null;
		try {
			SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			if(source==null)
				{
				in= srf.open(SamInputResource.of(System.in));
				}
			else
				{
				in= srf.open(SamInputResource.of(source));
				}
			iter = in.iterator();
			bindings.put("header",in.getFileHeader());
			bindings.put("iter",iter);
			bindings.put("format","sam");
			this.script.eval(bindings);
			return 0;
			} 
		catch (Exception e)
			{
			error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			CloserUtil.close(iter);
			bindings.remove("header");
			bindings.remove("iter");
			bindings.remove("format");
			}
		}
	
	public class Fasta 
		{
		private String sequence;
		private String name;
		public String getSequence() {
			return sequence;
			}
		public String getName() {
			return name;
			}
		public void print()
			{
			BioAlcidae.this.writer.print(">");
			BioAlcidae.this.writer.print(name);
			for(int i=0;i< sequence.length();++i)
				{
				if(i%60==0) BioAlcidae.this.writer.println();
				BioAlcidae.this.writer.print(this.charAt(i));
				}
			BioAlcidae.this.writer.println();
			}
		public int length()
			{
			return sequence.length();
			}
		public int getLength()
			{
			return sequence.length();
			}
		public int size()
			{
			return sequence.length();
			}
		public int getSize()
			{
			return sequence.length();
			}
		public char charAt(int i)
			{
			return sequence.charAt(i);
			}
		@Override
		public String toString() {
			return sequence;
			}
		}
	
	public class FastaIterator extends AbstractIterator<Fasta>
		{
		private LineIterator in=null;
		@Override
		protected Fasta advance()
			{
			Fasta f=null;
			for(;;)
				{
				if(!in.hasNext()) return null;
				String line=in.next();
				if(line.trim().isEmpty()) continue;
				if(!line.startsWith(">")) throw new RuntimeException("Expected '>'. Bad fasta line :"+line);
				f=new Fasta();
				f.name=line.substring(1);
				break;
				}
			StringBuilder sb=new StringBuilder();
			while(in.hasNext())
				{
				String line=in.peek();
				if(line.trim().isEmpty()) {in.next();continue;}
				if(line.startsWith(">")) break;//new sequence
				sb.append(in.next().trim());
				}
			f.sequence = sb.toString();
			return f;
			}
		}
	
	private int execute_fasta(String source) throws IOException
		{
		FastaIterator iter=new FastaIterator();
		try {
			iter.in = (source==null?IOUtils.openStdinForLineIterator():IOUtils.openURIForLineIterator(source));
			bindings.put("iter",iter);
			bindings.put("format","fasta");
			this.script.eval(bindings);
			return 0;
			} 
		catch (Exception e)
			{
			error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter.in);
			bindings.remove("iter");
			bindings.remove("format");
			}
		}
	
	private int execute_fastq(String source) throws IOException
		{
		InputStream in=null;
		FourLinesFastqReader iter=null;
		try {
			if(source==null)
				{
				in= System.in;
				}
			else
				{
				in= IOUtils.openURIForReading(source);
				}
			iter = new FourLinesFastqReader(in);
			bindings.put("iter",in);
			bindings.put("format","fastq");
			this.script.eval(bindings);
			return 0;
			} 
		catch (Exception e)
			{
			error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(in);
			bindings.remove("header");
			bindings.remove("iter");
			bindings.remove("format");
			}
		}
	
	
	private int execute(FORMAT fmt,String source) throws IOException
		{
		switch(fmt)
			{
			case VCF: return execute_vcf(source);
			case BAM: case SAM: return execute_bam(source);
			case FASTQ: return execute_fastq(source);
			case FASTA: return execute_fasta(source);
			default: throw new IllegalStateException();
			}
		}
	
	@Override
	public int executeKnime(List<String> args)
		{
		try
			{
			this.bindings = this.engine.createBindings();
			if(getOutputFile()==null)
				{
				this.writer = new PrintWriter(System.out);
				}
			else
				{
				this.writer = new PrintWriter(getOutputFile());
				}
			this.bindings.put("out",this.writer);
			
			if(args.isEmpty()  )
				{
				if(this.format==null)
					{
					error("Format must be specified when reading from stdin");
					return -1;
					}
				return execute(this.format,null);
				}
			for(String filename:IOUtils.unrollFiles(args))
				{
				this.bindings.put("FILENAME",filename);
				FORMAT fmt=this.format;
				if(fmt==null)
					{
					for(FORMAT t:FORMAT.values())
						{
						if(t.canAs(filename))
							{
							fmt=t;
							break;
							}
						}
					}
				if(fmt==null)
					{
					error("Cannot get file format for "+filename);
					return -1;
					}
				if( execute(fmt,filename) !=0 )
					{
					error("Failure for "+filename);
					return -1;
					}
				}
				
			
			System.out.flush();
			return 0;
			}
		catch(Exception err)
			{
			error(err);	
			return -1;
			}
		finally
			{
			if(this.writer!=null) this.writer.flush();
			if(getOutputFile()!=null)
				{
				CloserUtil.close(this.writer);
				}
			this.writer=null;
			bindings=null;
			}
		}
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -e (js expression). Optional.");
		out.println(" -f (js file). Optional.");
		out.println(" -F (format) force format: one of "+Arrays.asList(FORMAT.values()));
		super.printOptions(out);
		}
	
	
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:e:f:F:"))!=-1)
			{
			switch(c)
				{
				case 'F':this.setFormat(opt.getOptArg());break;
				case 'o':this.setOutputFile(new File(opt.getOptArg()));break;
				case 'e':this.setScriptExpression(opt.getOptArg());break;
				case 'f':this.setScriptFile(new File(opt.getOptArg()));break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
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
	
	public static void main(String[] args) {
		new BioAlcidae().instanceMainWithExit(args);
	}
}
