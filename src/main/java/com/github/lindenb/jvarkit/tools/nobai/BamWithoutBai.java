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
package com.github.lindenb.jvarkit.tools.nobai;

import java.io.ByteArrayInputStream;
import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;
import java.io.SequenceInputStream;
import java.nio.file.Path;
import java.util.List;

import org.apache.http.impl.client.HttpClientBuilder;
import org.apache.http.impl.conn.DefaultProxyRoutePlanner;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;

import org.apache.http.HttpEntity;
import org.apache.http.HttpHost;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;

/**
BEGIN_DOC

Motivation: query a remote bam without bai.

END_DOC

**/

@Program(name="bamwithoutbai",
description="Query a Remote BAM without bai",
keywords={"bam","sam","bai"},
creationDate="20191213",
modificationDate="20191213",
generate_doc=false
)
public class BamWithoutBai extends Launcher{
	private static final Logger LOG = Logger.build(BamWithoutBai.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-r","--region","--interval"},description=IntervalParserFactory.OPT_DESC,required=true)
	private String intervalStr ="";
	@Parameter(names={"--buffer-size"},description="Buffer size. Must be large enough to contains a few reads and a BGZF block")
	private int byte_buffer_size=100_000;
	@Parameter(names={"--repeat"},description="Dichotomy repeat")
	private int dichotomy_repeat=20;

	
	private CloseableHttpClient httpClient = null;
	
	/** empty closeable iterator */
	private static final CloseableIterator<SAMRecord> EOF_ITER  = new CloseableIterator<SAMRecord>() {
		@Override
		public SAMRecord next() { return null; }
		@Override
		public boolean hasNext() { return false; }
		@Override
		public void close() { }
	};
	

		
	
	private class MySamRecordIterator 
		extends AbstractIterator<SAMRecord>
		implements CloseableIterator<SAMRecord> {
        private final BAMRecordCodec bamRecordCodec;
        private Runnable onClose = null;
        MySamRecordIterator(final SAMFileHeader header,final InputStream in,final String url) {
        	this.bamRecordCodec = new BAMRecordCodec(header);
        	this.bamRecordCodec.setInputStream(in,url);
        }
        
		@Override
		protected SAMRecord advance() {
            final SAMRecord next = this.bamRecordCodec.decode();
            if(next==null) close();
            return next;
			}
		@Override
		public void close() {
			if(onClose!=null) onClose.run();
			}
		}
	
	private class MySamFilterIterator 
		extends AbstractIterator<SAMRecord>
		implements CloseableIterator<SAMRecord> {
		final CloseableIterator<SAMRecord> delegate;
		final QueryInterval interval;
		MySamFilterIterator(final CloseableIterator<SAMRecord> delegate,final QueryInterval interval) {
	    	this.delegate = delegate;
			this.interval = interval;
	    	
	    }
	    
		@Override
		protected SAMRecord advance() {
			while(delegate.hasNext()) {
				final SAMRecord rec = delegate.next();
				if(rec.getReferenceIndex()!=this.interval.referenceIndex) continue;
				if(!CoordMath.overlaps(rec.getStart(), rec.getEnd(), this.interval.start, this.interval.end)) continue;
				return rec;
				}
			return null;
			}
		@Override
		public void close() {
			delegate.close();
			}
		}

	
	
	private class CustomSamReader implements Closeable {
		final String url;
		final long content_length;
		final HttpGet httpGet;
		final SAMFileHeader samFileHeader;

		
		CustomSamReader(final String url ) throws IOException {
				this.url=url;
				this.httpGet = new HttpGet(url);
				
				final SamReaderFactory srf = SamReaderFactory.makeDefault().
						validationStringency(ValidationStringency.LENIENT);
				
				// get Header
				try( CloseableHttpResponse httpResponse = httpClient.execute(httpGet)) {
					 final int responseCode = httpResponse.getStatusLine().getStatusCode();
					  
					 if(responseCode != 200)
					 	{
						throw new RuntimeIOException("Response code was not 200. Detected response was "+responseCode);
					 	}
					 final HttpEntity entity = httpResponse.getEntity();
					 this.content_length = entity.getContentLength();
					 try(InputStream in=entity.getContent()) {
						try(SamReader sr= srf.open(SamInputResource.of(in))) {
							this.samFileHeader = sr.getFileHeader();
						}
					 }
				}
				if(!this.samFileHeader.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
					throw new IOException("Bam "+url+" is not sorted on coordinate: "+this.samFileHeader.getSortOrder());
				}
			}
		
		public SAMFileHeader getFileHeader() {
			return samFileHeader;
		}

		
		private CloseableIterator<SAMRecord> queryAtOffset(final long start_offset) throws IOException {
			if(start_offset>=this.content_length) return EOF_ITER;
			LOG.warn("opening at offset="+start_offset);
			final byte buffer[]  =  new byte[byte_buffer_size];
            httpGet.setHeader("Range", "bytes="+start_offset+"-"+this.content_length);
            CloseableHttpResponse httpResponse = httpClient.execute(httpGet);
            final int responseCode = httpResponse.getStatusLine().getStatusCode();

            if(responseCode != 200)
 			 	{
            	httpResponse.close();
            	throw new RuntimeIOException("Response code was not 200 for offset '"+start_offset+"'. Detected response was "+responseCode);
 			 	}
            	
			 final HttpEntity entity = httpResponse.getEntity();
			 final InputStream in;
			 try {
				in = entity.getContent();
			 	}
			 catch(final Throwable err) {
				httpResponse.close();
				throw new RuntimeIOException(err);
			 	}
			 //fill buffer
			 int nBytesRead = 0;
			 while(nBytesRead< buffer.length) {
				int n=in.read(buffer,nBytesRead,buffer.length-nBytesRead);
				if(n==-1) break;
				nBytesRead+=n;
			 	}
			 if(nBytesRead==0) {
				 LOG.warn("nothing read at "+start_offset);
				 in.close();
				 return EOF_ITER;
			 }
				 
				 
			 MySamRecordIterator mySamRecordIter = null;
			 for(int offset=0;offset < buffer.length;offset++) {
				 mySamRecordIter = null;
				 final ByteArrayInputStream bais = new ByteArrayInputStream(buffer, 0, buffer.length);
				 final SequenceInputStream mergeIn= new SequenceInputStream(bais,in);
				 final BlockCompressedInputStream bci = new BlockCompressedInputStream(mergeIn);
				 try {
					mySamRecordIter = new MySamRecordIterator(this.samFileHeader,bci,this.url);
					mySamRecordIter.hasNext();//force checking next
					mySamRecordIter.onClose=()->{
						try {
							in.close();
							httpResponse.close();
							}
						catch(final IOException err) 
							{
							LOG.warning(err);
							}
						};
					return mySamRecordIter;
				 	}
				 catch(final Throwable err) {
					 System.err.println("offset:"+offset+" "+err.getMessage());
					 mySamRecordIter = null;
				 	}
			 	}
			in.close();
			httpResponse.close();
			return EOF_ITER;
			}
		
		
		
		public CloseableIterator<SAMRecord> query(final Locatable locatable) throws IOException {
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(samFileHeader);
			final int tid = dict.getSequenceIndex(locatable.getContig());
			if(tid==-1) {
				LOG.warn("unknown contig "+locatable.getContig());
				return EOF_ITER;
			}
			final QueryInterval userQueryInterval = new QueryInterval(tid, locatable.getStart(),locatable.getEnd());
			long byte_start = 0L;
			long byte_end=this.content_length;
			int repeat_dichotomy = dichotomy_repeat;
		    long len = byte_start - byte_end;
		    while (len > 0L)
		            {
		    		repeat_dichotomy--;
		            final long half = len / 2;
		            final long middle = byte_start + half;
		            CloseableIterator<SAMRecord> iter1 = queryAtOffset(middle);
		            if(iter1.hasNext())
		            	{
		            	if(repeat_dichotomy==0) return new MySamFilterIterator(iter1,userQueryInterval);
		            	final SAMRecord rec = iter1.next();
		            	
		            	LOG.debug("got "+rec.getSAMString());
		            	if(
		            		(rec.getReferenceIndex()!=SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && rec.getReferenceIndex()< userQueryInterval.referenceIndex) ||
		            		(rec.getReferenceIndex()== userQueryInterval.referenceIndex && rec.getAlignmentStart() <= userQueryInterval.start)
		            		)
		            		{
		            		byte_start = middle + 1;
		                    len = len - half - 1;
		            		}
		            	else
		            		{
		            		len = half;
		            		}
		            	
		            	}
		            else
		            	{
		            	return EOF_ITER;
		            	}
		            iter1.close();
		            }
				
			return EOF_ITER;
			}
		
		
		@Override
		public void close()  {
			
			}
	}
	
	
	
	@Override
	public int doWork(final List<String> args) {
	
	try {
		final String urlStr = oneAndOnlyOneFile(args);
		if(!( urlStr.startsWith("http://") || urlStr.startsWith("https://"))) {
			LOG.error("No a remote http url: "+urlStr);
			return -1;
		}


		HttpHost proxy = new HttpHost("cache.ha.univ-nantes.fr", 3128, "http");
		/** create http client */
		this.httpClient = HttpClientBuilder.
				create().
				setRoutePlanner(new DefaultProxyRoutePlanner(proxy)).
				setUserAgent(IOUtils.getDefaultUserAgent()).
				build();
		
		
		LOG.debug("getting header for "+urlStr);
		try( CustomSamReader sr =  new CustomSamReader(urlStr)) {
			final SAMFileHeader samFileHeader = sr.getFileHeader();
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(samFileHeader);
			final Locatable userInterval = IntervalParserFactory.newInstance(dict).make().apply(this.intervalStr).orElse(null);
			if(userInterval==null) {
				LOG.error("cannot parse interval "+this.intervalStr+" for "+urlStr);
				return -1;
				}
			try(CloseableIterator<SAMRecord> iter = sr.query(userInterval)) {
				while(iter.hasNext()) {
					final SAMRecord rec=iter.next();
					System.err.println(""+rec);
					}
				}
			}
		
		return 0;
	} catch(final Throwable err) {
		LOG.error(err);
		return -1;
	} finally {
		CloserUtil.close(this.httpClient);
	}
	
	}
	
public static void main(final String[] args) {
	new BamWithoutBai().instanceMainWithExit(args);
	}
}
