/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.ContigPosRef;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

public abstract class AbstractVCFCompareBase extends Launcher
	{
	private final Logger LOG=Logger.build(AbstractVCFCompareBase.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	protected File outputFile = null;

	@ParametersDelegate
	protected WritingSortingCollection sortingCollectionArgs= new WritingSortingCollection(); 

	/** input files */
	protected List<Input> inputs=new ArrayList<Input>();

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
		public ContigPosRef getContigPosRef() 
			{
			return new ContigPosRef(getContig(), getStart(), getReference());
			}
		public String getChrom()
			{	
			return getContig();
			}
		public String getContig()
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
		public LineAndFile decode(final DataInputStream dis) throws IOException
			{
			final LineAndFile v=new LineAndFile();
			try {
				v.fileIdx=dis.readInt();
			} catch (Exception e) {
				return null;
				}
			v.line= readString(dis);
			return v;
			}
		@Override
		public void encode(final DataOutputStream dos, final LineAndFile v)
				throws IOException
			{
			dos.writeInt(v.fileIdx);
			writeString(dos,v.line);
			}
		@Override
		public LineAndFileCodec clone() {
			return new LineAndFileCodec();
			}
		}
	
	
	
	protected static class LineAndFileComparator implements Comparator<LineAndFile>
		{
		@Override
		public int compare(final LineAndFile v1,final LineAndFile v2)
			{
			int i= v1.getContigPosRef().compareTo(v2.getContigPosRef());
			if(i!=0) return i;
			i =  v1.fileIdx - v2.fileIdx;
			if(i!=0) return i;
			return v1.line.compareTo(v2.line);
			}
		}
	
	
	
	protected AbstractVCFCompareBase()
		{
		}
	
/*	
	
	private Comparator<LineAndFile> createLineAndFileComparator()
		{
		return new LineAndFileComparator();
		}
*/
	
	protected Input createInput()
		{
		return new Input();
		}
	
	/** creates a new input=(filename, codec , vcfheader ) */
	private Input createInput(final String filename,final List<String> headerLines)
		{
		final Input input=createInput();
		input.filename=filename;
		input.codecAndHeader = VCFUtils.parseHeader(headerLines);
		return input;
		}
	
	/** give a chance to simplify the vcf line. returned null string will be ignored */
	protected String simplify(String line,int input_idx)
		{
		return line;
		}
	
	/** insert all  variant of vcfUri into the sorting collection */
	protected Input put(final SortingCollection<LineAndFile> variants, String vcfUri)
		throws IOException
		{
		LOG.info("begin inserting "+vcfUri);
		LineIterator iter=null;
		if(vcfUri==null)
			{
			vcfUri="stdin";
			iter=IOUtils.openStreamForLineIterator(stdin());
			}
		else
			{
			iter=IOUtils.openFileForLineIterator(new File(vcfUri));
			}
		
		final List<String> headerLines=new ArrayList<String>();
		while(iter.hasNext() && iter.peek().startsWith("#"))
			{
			headerLines.add(iter.next());
			}
		if(headerLines.isEmpty()) throw new IOException("Not header found in "+vcfUri);
		final Input input=createInput(vcfUri, headerLines);
		input.file_id=this.inputs.size();
		this.inputs.add(input);
		final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(input.codecAndHeader.header.getSequenceDictionary());
		while(iter.hasNext())
			{
			final String line=iter.next();
			final LineAndFile laf=new LineAndFile();
			laf.fileIdx=input.file_id;
			laf.line=simplify(line,laf.fileIdx);
			if(laf.line==null) continue;
			progress.watch(laf.getChrom(), laf.getStart());
			variants.add(laf);
			input.count++;
			}
		progress.finish();
		LOG.info("end inserting "+vcfUri+" N="+input.count);
		CloserUtil.close(iter);
		return input;
		}
	}
