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

*/
package com.github.lindenb.jvarkit.tools.vcfviewgui;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.Optional;

import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;

public class BamFile implements NgsFile<SAMFileHeader,SAMRecord>{
    private static final Logger LOG= Logger.build(BamFile.class).make();

	private final SamReader samReader;
	private final String source;
	private final Optional<File> indexFile;
	private boolean delete_index_on_close=true;
	private BamFile(final String source, final SamReader samReader,final Optional<File> indexFile)
		throws IOException
		{
		this.samReader = samReader;
		this.source = source;
		this.indexFile = indexFile;
		
	    if(!this.samReader.hasIndex())
       	{
       	this.samReader.close();
       	throw new IOException("Bam without index "+source);
       	}
       
       if(this.samReader.getFileHeader()==null)
   		{
	    	this.samReader.close();
	    	throw new IOException("Bam without header "+source);
	    	}
       if(this.samReader.getFileHeader().getSequenceDictionary()==null)
			{
	    	this.samReader.close();
	    	throw new IOException("Bam without dictionary "+source);
	    	}

		}
	
	public static BamFile newInstance(final String s) throws IOException {
		return IOUtil.isUrl(s)?
				newInstance(new URL(s)):
				newInstance(new File(s))
				;
		}
	
	public static BamFile newInstance(final File f) throws IOException {
		IOUtil.assertFileIsReadable(f);
		final SamReaderFactory srf = SamReaderFactory.makeDefault();
		srf.validationStringency(ValidationStringency.LENIENT);
		return new BamFile(f.getAbsolutePath(),srf.open(f),Optional.empty());
		}
	
	public static BamFile newInstance(final URL url) throws IOException {		
    	final File baiFile=File.createTempFile("tmp.", BAMIndex.BAMIndexSuffix);
    	Optional<File> savedBaiFile=Optional.empty();
    	baiFile.deleteOnExit();
    	for(int i=0;i< 2;++i)
    		{
    		final String bamurl=url.toExternalForm();
    		final String baiurl=(i==0?
    				bamurl:bamurl.substring(0, bamurl.length()-BamFileIoUtils.BAM_FILE_EXTENSION.length())
    				)+ BAMIndex.BAMIndexSuffix;
    		InputStream in=null;
    		FileOutputStream out=null;
    		try {
    			LOG.info("trying "+baiurl);
				in = new URL(baiurl).openStream();
				out = new FileOutputStream(baiFile);
				IOUtil.copyStream(in, out);
				out.flush();
				out.close();
				in.close();
				savedBaiFile = Optional.of(baiFile);
				break;
			} catch (final IOException err) {
				baiFile.delete();
				LOG.info("Cannot fetch "+baiurl+" : "+err.getMessage());
				}
    		finally
    			{
    			CloserUtil.close(out);
    			CloserUtil.close(in);
    			}
    		}
    	if(!savedBaiFile.isPresent())
    		{
    		throw new IOException("cannot get a bam index file for "+url);
    		}	
		final SamReaderFactory srf = SamReaderFactory.makeDefault();
		srf.validationStringency(ValidationStringency.LENIENT);
		final SamInputResource sir= SamInputResource.of(url);
		sir.index(savedBaiFile.get());
		return new BamFile(url.toExternalForm(),srf.open(sir),savedBaiFile);
		}
	
	@Override
	public SAMFileHeader getHeader() {
		return samReader.getFileHeader();
	}

	@Override
	public SAMSequenceDictionary getSequenceDictionary() {
		return getHeader().getSequenceDictionary();
	}

	
	public CloseableIterator<SAMRecord> queryUnmapped()  throws IOException {
		return samReader.queryUnmapped();
	}
	
	@Override
	public CloseableIterator<SAMRecord> iterator() throws IOException {
		return samReader.iterator();
	}

	@Override
	public CloseableIterator<SAMRecord> iterator(String contig, int start, int end) throws IOException {
		return samReader.query(contig,start,end,false);
	}

	@Override
	public void close() {
		if(this.indexFile.isPresent() && this.delete_index_on_close)
			{
			this.indexFile.get().delete();	
			}
		CloserUtil.close(this.samReader);
		}

	@Override
	public String getSource() {
		return this.source;
	}

	@Override
	public BamFile reOpen() throws IOException {
		final String url=this.getSource();
		final SamReaderFactory srf = SamReaderFactory.makeDefault();
		srf.validationStringency(ValidationStringency.LENIENT);
		final SamInputResource sir;
		if(IOUtil.isUrl(url))
			{
			sir = SamInputResource.of(new URL(url));
			if(!this.indexFile.isPresent()) throw new IOException("Boum");
			sir.index(this.indexFile.get());
			}
		else
			{
			sir = SamInputResource.of(new File(url));
			}
		
		final BamFile bf = new BamFile(url,srf.open(sir),this.indexFile);
		bf.delete_index_on_close=false;
		return bf;
		}
}
