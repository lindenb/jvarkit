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

import java.io.Closeable;
import java.io.File;
import java.io.IOException;

import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;



public class TribbleVcfFileReader
	implements Closeable
	{
	private static final Logger LOG=Logger.build(TribbleVcfFileReader.class).make();
	private FeatureCodec<VariantContext, LineIterator> tribbleCodec=null;
	private Index tribbleIndex=null;
    private File source;
    private AbstractFeatureReader<VariantContext,LineIterator> reader;
    
    
    public TribbleVcfFileReader(final File vcf) throws IOException
    	{
    	this.source=vcf;
    	if(vcf==null) throw new NullPointerException("vcf file==null");
    	IOUtil.assertFileIsReadable(vcf);
    	this.tribbleCodec = VCFUtils.createAsciiFeatureCodec();
    	
    	File indexFile=Tribble.indexFile(this.source);
    	
    	if(indexFile.exists() && indexFile.lastModified()> vcf.lastModified())
		 	{
			LOG.info("loading index in memory for "+this.source+" index="+indexFile);
			this.tribbleIndex=IndexFactory.loadIndex(indexFile.getPath());
		 	}
    	else
		 	{
			LOG.info("create index from file "+this.source);
			this.tribbleIndex=IndexFactory.createLinearIndex(vcf, tribbleCodec);
		 	}
    	this.reader =
    			AbstractFeatureReader.getFeatureReader(
    					vcf.getPath(),
    					this.tribbleCodec,
    					this.tribbleIndex
    					);
    	
    	}
    
    public File getSource()
    	{
    	return this.source;
    	}
    
	public VCFHeader getHeader()
		{
	    return VCFHeader.class.cast(this.reader.getHeader());
		}
    
    public CloseableTribbleIterator<VariantContext> iterator(String chrom,int start,int end)
    	throws IOException
    	{
    	if(this.reader==null) throw new IllegalStateException("Reader==null");
    	return this.reader.query(chrom, start, end);
    	}
    
    public CloseableTribbleIterator<VariantContext> iterator()
        	throws IOException
    	{
    	if(this.reader==null) throw new IllegalStateException("Reader==null");
    	return this.reader.iterator();
    	}
    @Override
    public void close() throws IOException {
    	CloserUtil.close(this.reader);
    	this.reader=null;
    	this.tribbleCodec=null;
    	}
    
    
	}
