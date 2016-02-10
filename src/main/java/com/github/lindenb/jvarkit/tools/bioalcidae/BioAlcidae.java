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

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.List;

import javax.script.Bindings;
import javax.script.CompiledScript;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/**
 * BioAlcidae
 *
 */
public class BioAlcidae
	extends AbstractBioAlcidae
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(BioAlcidae.class);
	
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
		private FORMAT format= null;
		private PrintWriter writer=null;


		

	
	/** moved to public for knime */
	public  Collection<Throwable> executeAsVcf(String source) throws IOException
		{
		LOG.info("source: "+source);
		VcfIterator in=null;
		try {
			
			return executeAsVcf(source);
			} 
		catch (Exception e)
			{
			return wrapException(e);
			}
		finally
			{
			CloserUtil.close(in);
			}
		}
	
	public  Collection<Throwable> executeAsVcf(final VcfIterator in) throws IOException
		{
		try {
			bindings.put("codec",in.getCodec());
			bindings.put("header",in.getHeader());
			bindings.put("iter",in);
			bindings.put("format","vcf");
			this.script.eval(bindings);
			this.writer.flush();
			return RETURN_OK;
			} 
		catch (Exception e)
			{
			return wrapException(e);
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

	
	
	private  Collection<Throwable> execute_bam(String source) throws IOException
		{
		SamReader in=null;
		SAMRecordIterator iter=null;
		try {
			SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			if(source==null)
				{
				in= srf.open(SamInputResource.of(stdin()));
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
			this.writer.flush();
			return RETURN_OK;
			} 
		catch (Exception e)
			{
			LOG.error(e);
			return wrapException(e);
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
	
	private  Collection<Throwable> execute_fasta(String source) throws IOException
		{
		FastaIterator iter=new FastaIterator();
		try {
			iter.in = (source==null?IOUtils.openStdinForLineIterator():IOUtils.openURIForLineIterator(source));
			bindings.put("iter",iter);
			bindings.put("format","fasta");
			this.script.eval(bindings);
			this.writer.flush();
			return RETURN_OK;
			} 
		catch (Exception e)
			{
			LOG.error(e);
			return wrapException(e);
			}
		finally
			{
			CloserUtil.close(iter.in);
			bindings.remove("iter");
			bindings.remove("format");
			}
		}
	
	private  Collection<Throwable> execute_fastq(String source) throws IOException
		{
		InputStream in=null;
		FourLinesFastqReader iter=null;
		try {
			if(source==null)
				{
				in= stdin();
				}
			else
				{
				in= IOUtils.openURIForReading(source);
				}
			iter = new FourLinesFastqReader(in);
			bindings.put("iter",in);
			bindings.put("format","fastq");
			this.script.eval(bindings);
			this.writer.flush();
			return RETURN_OK;
			} 
		catch (Exception e)
			{
			LOG.error(e);
			return wrapException(e);
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
	
	
	private  Collection<Throwable> execute(FORMAT fmt,String source) throws IOException
		{
		switch(fmt)
			{
			case VCF: return executeAsVcf(source);
			case BAM: case SAM: return execute_bam(source);
			case FASTQ: return execute_fastq(source);
			case FASTA: return execute_fasta(source);
			default: return wrapException(new IllegalStateException());
			}
		}
	
	@Override
	public Collection<Throwable> initializeKnime()
		{
		if(this.formatString!=null)
			{
			try {
				this.format=FORMAT.valueOf(this.formatString.toUpperCase());
				} catch (Exception e) {
				return wrapException(e);
				}
			}
		try
			{
			this.script = super.compileJavascript();
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		return super.initializeKnime();
		}
	
	@Override
	public void disposeKnime() {
		CloserUtil.close(this.writer);
		this.writer=null;
		this.bindings=null;
		this.bindings=null;
		this.format=null;
		super.disposeKnime();
		}
	
	/** moved to public for knime */
	public void initializeJavaScript() throws IOException
		{
		this.bindings = this.script.getEngine().createBindings();
		this.writer = super.openFileOrStdoutAsPrintWriter();
		this.bindings.put("out",this.writer);
		}
	
	@Override
	public Collection<Throwable> call() throws Exception
		{
		final List<String> args = getInputFiles();
		try
			{
			this.initializeJavaScript();
			
			if(args.isEmpty()  )
				{
				if(this.format==null)
					{
					return wrapException("Format must be specified when reading from stdin");
					}
				return execute(this.format,null);
				}
			for(final String filename:IOUtils.unrollFiles(args))
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
					return wrapException("Cannot get file format for "+filename);
					}
				Collection<Throwable> errors= execute(fmt,filename);
				if(!errors.isEmpty())
					{
					return errors;
					}
				}
				
			
			this.writer.flush();
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			if(this.writer!=null) this.writer.flush();
			if(getOutputFile()!=null)
				{
				CloserUtil.close(this.writer);
				}
			this.writer=null;
			this.bindings=null;
			this.bindings=null;
			this.format=null;
			}
		}
	
	
	
	
	
	
	
	public static void main(String[] args) {
		new BioAlcidae().instanceMainWithExit(args);
	}
}
