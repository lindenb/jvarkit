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
* 2015 creation

*/

package com.github.lindenb.jvarkit.util.vcf;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/* Tabix or Tribble indexed file reader */
public class IndexedVcfFileReader
	implements Closeable
	{
	private enum Type {Tabix,Tribble};
	private Closeable reader;
	private Type type;
	private Object source;
	public IndexedVcfFileReader(File vcf) throws IOException
		{
		init(vcf);
		}
	
	private void init(File vcf)  throws IOException
		{
    	if(vcf==null) throw new NullPointerException("vcf file==null");
    	if(!vcf.isFile())  throw new IOException("vcf is not a file "+vcf);
    	if(!vcf.canRead())  throw new IOException("cannot read "+vcf);
    	if(vcf.getName().endsWith(".gz"))
    		{
    		this.type=Type.Tabix;
    		this.reader=new TabixVcfFileReader(vcf.getPath());
    		}
    	else
    		{
    		this.type=Type.Tribble;
    		this.reader=new TribbleVcfFileReader(vcf);
    		}
    	this.source=vcf;
		}
	
	public IndexedVcfFileReader(String onlyTabixOrLocal) throws IOException
		{
		if(onlyTabixOrLocal==null) throw new NullPointerException("vcf file==null");
		
		if(IOUtil.isUrl(onlyTabixOrLocal))
			{
			this.type=Type.Tabix;
			this.reader=new TabixVcfFileReader(onlyTabixOrLocal);
			this.source=onlyTabixOrLocal;
			}
		else
			{
			init(new File(onlyTabixOrLocal));
			}
		}

	
	private void checkOpen()
		{
		if(this.reader==null)
				throw new IllegalStateException("vcf reader is closed "+getSource());
		}
	
	/* File or String */
	public Object getSource()
		{
		return this.source;
		}
	
	public VCFHeader getHeader()
		{
		checkOpen();
		switch(type)
			{
			case Tabix: return TabixVcfFileReader.class.cast(reader).getHeader();
			case Tribble: return TribbleVcfFileReader.class.cast(reader).getHeader();
			default: throw new IllegalStateException();
			}
		}
	public CloseableIterator<VariantContext>
		iterator(String chrom,int start,int end)
		throws IOException
		{
		checkOpen();
		switch(type)
			{
			case Tabix: return new CIter(TabixVcfFileReader.class.cast(reader).iterator(chrom, start, end));
			case Tribble: return TribbleVcfFileReader.class.cast(reader).iterator(chrom,start,end);
			default: throw new IllegalStateException();
			}
		}
	
	public List<VariantContext> getVariants(String chrom,int start,int end) throws IOException
		{
		checkOpen();
		CloseableIterator<VariantContext> iter=null;
		try
			{
			List<VariantContext> L= new ArrayList<>();
			iter = iterator(chrom,start,end);
			while(iter.hasNext())
				{
				L.add(iter.next());
				}
			return L;
			}
		finally
			{
			CloserUtil.close(iter);
			}
		}
	
	
	private static class CIter extends AbstractIterator<VariantContext>
		implements CloseableIterator<VariantContext>
		{
		Iterator<VariantContext> delegate=null;
		CIter(Iterator<VariantContext> delegate)
			{
			this.delegate=delegate;
			}
		@Override
		protected VariantContext advance() {
			return this.delegate!=null && this.delegate.hasNext()?this.delegate.next():null;
			}
		@Override
		public void close() {
			CloserUtil.close(this.delegate);
			this.delegate=null;
			}
		}
	
	@Override
	public void close() throws IOException {
		CloserUtil.close(reader);
		reader=null;
		}
	}
