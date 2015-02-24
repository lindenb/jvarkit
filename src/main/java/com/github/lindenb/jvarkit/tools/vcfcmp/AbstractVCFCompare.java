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
package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReaderUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.cli.GetOpt;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SortingCollectionFactory;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

public abstract class AbstractVCFCompare extends AbstractCommandLineProgram
	{
	/** input files */
	protected List<Input> inputs=new ArrayList<Input>();
	protected SortingCollectionFactory<LineAndFile> factory=new SortingCollectionFactory<LineAndFile>();

	/** filename, codec, header, file-id, count */
	protected class Input
		{
		String filename;
		VCFUtils.CodecAndHeader codecAndHeader;
		int file_id=-1;
		long count=0L;
		}
	
	/** line associated to file */
	protected class LineAndFile
		{
		int fileIdx;
		String line;
		
		private VariantContext _ctx=null;
		VariantContext getContext()
			{
			if(this._ctx==null)
				{
				this._ctx=inputs.get(this.fileIdx).codecAndHeader.codec.decode(this.line);
				}
			return this._ctx;
			}
		private String token(int index)
			{
			int i=0;
			int prev=0;
			while(prev< line.length())
				{
				int tab=this.line.indexOf('\t',prev);
				if(tab==-1)
					{
					tab=this.line.length();
					}
				if(i==index) return this.line.substring(prev, tab);
				prev=tab+1;
				++i;
				}
			throw new IllegalStateException("Bad VCF line :"+line+" cannot get tokens["+index+"]");
			}
		
		public String getChrom()
			{	
			return token(0);
			}
		public int getStart()
			{	
			return Integer.parseInt(token(1));
			}
		public Allele getReference()
			{	
			return Allele.create(token(3),true);
			}
		}
	
	protected class LineAndFileCodec extends AbstractDataCodec<LineAndFile>
		{
		@Override
		public LineAndFile decode(DataInputStream dis) throws IOException
			{
			LineAndFile v=new LineAndFile();
			try {
				v.fileIdx=dis.readInt();
			} catch (Exception e) {
				return null;
				}
			v.line= readString(dis);
			return v;
			}
		@Override
		public void encode(DataOutputStream dos, LineAndFile v)
				throws IOException
			{
			dos.writeInt(v.fileIdx);
			writeString(dos,v.line);
			}
		@Override
		public AbstractDataCodec<LineAndFile> clone() {
			return new LineAndFileCodec();
			}
		}

	
	
	
	
	
	protected class LineAndFileComparator implements Comparator<LineAndFile>
		{
		@Override
		public int compare(LineAndFile v1, LineAndFile v2)
			{
			int i=v1.getChrom().compareTo(v2.getChrom());
			if(i!=0) return i;
			i=v1.getStart()-v2.getStart();
			if(i!=0) return i;
			i=v1.getReference().compareTo(v2.getReference());
			if(i!=0) return i;
			return 0;
			}
		}
	
	
	
	protected AbstractVCFCompare()
		{
		}
	
	
	@Override
	public String getProgramDescription()
		{
		return "Compares two VCF files.";
		}
	
	@Override
	public void printOptions(PrintStream out) {
		out.println(" -M (int) Max records in RAM. Optional.");
		out.println(" -T (dir) add temporary directory. Optional");
		super.printOptions(out);
		}
	
	@Override
	protected String getGetOptDefault() {
		return super.getGetOptDefault()+"M:T:";
		}
	
	protected Comparator<LineAndFile> createLineAndFileComparator()
		{
		return new LineAndFileComparator();
		}

	
	protected Input createInput()
		{
		return new Input();
		}
	
	private Input createInput(String filename,List<String> headerLines)
		{
		Input input=createInput();
		input.filename=filename;
		input.codecAndHeader = VCFUtils.parseHeader(headerLines);
		return input;
		}
	
	protected String simplify(String line,int input_idx)
		{
		return line;
		}
	
	protected Input put(SortingCollection<LineAndFile> variants, String vcfUri)
		throws IOException
		{
		info("begin inserting "+vcfUri);
		LineIterator iter=null;
		if(vcfUri==null)
			{
			vcfUri="stdin";
			iter=new LineIteratorImpl(LineReaderUtil.fromBufferedStream(System.in));
			}
		else
			{
			iter=IOUtils.openFileForLineIterator(new File(vcfUri));
			}
		
		List<String> headerLines=new ArrayList<String>();
		while(iter.hasNext() && iter.peek().startsWith("#"))
			{
			headerLines.add(iter.next());
			}
		if(headerLines.isEmpty()) throw new IOException("Not header found in "+vcfUri);
		Input input=createInput(vcfUri, headerLines);
		input.file_id=this.inputs.size();
		this.inputs.add(input);
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(input.codecAndHeader.header.getSequenceDictionary());
		while(iter.hasNext())
			{
			String line=iter.next();
			LineAndFile laf=new LineAndFile();
			laf.fileIdx=input.file_id;
			laf.line=simplify(line,laf.fileIdx);
			progress.watch(laf.getChrom(), laf.getStart());
			variants.add(laf);
			input.count++;
			}
		progress.finish();
		info("end inserting "+vcfUri+" N="+input.count);
		CloserUtil.close(iter);
		return input;
		}
	
	@Override
	protected GetOptStatus
		handleOtherOptions(int c, GetOpt opt, String[] args)
		{
		switch(c)
			{
			case 'M': factory.setMaxRecordsInRAM(Math.max(1,Integer.parseInt(opt.getOptArg()))); return GetOptStatus.OK;
			case 'T': super.addTmpDirectory(new File(opt.getOptArg())) ; return GetOptStatus.OK;
			default: return super.handleOtherOptions(c, opt, args);
			}
		}
	
	}
