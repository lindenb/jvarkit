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

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;
import java.io.SequenceInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
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
import htsjdk.samtools.seekablestream.SeekableBufferedStream;
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

## Example

'https://www.encodeproject.org/files/ENCFF741DEO/@@download/ENCFF741DEO.bam'

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
	@Parameter(names={"--debug"},description="Enable debugging.")
	private boolean do_debug = false;

	
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
	

		
	private class HeadingSamRecordIterator extends  AbstractIterator<SAMRecord>
		implements CloseableIterator<SAMRecord> {
        private SamReader reader;
        private CloseableIterator<SAMRecord> delegate;
        private CloseableHttpResponse httpResponse;
        HeadingSamRecordIterator(final HttpGet httpGet) throws IOException{
            this.httpResponse = httpClient.execute(httpGet);
            final int responseCode = httpResponse.getStatusLine().getStatusCode();

            if(responseCode != 200)
 			 	{
            	this.httpResponse.close();
            	throw new RuntimeIOException("Response code was not 200 . Detected response was "+responseCode);
 			 	}
            	
			 final HttpEntity entity = httpResponse.getEntity();
			final SamReaderFactory srf = SamReaderFactory.makeDefault().
					validationStringency(ValidationStringency.LENIENT);
			this.reader = srf.open(SamInputResource.of(entity.getContent()));
			this.delegate = this.reader.iterator();
        }
        
		@Override
		protected SAMRecord advance() {
            return delegate.hasNext()?delegate.next():null;
			}
		@Override
		public void close() {
			CloserUtil.close(delegate);
			CloserUtil.close(reader);
			CloserUtil.close(httpResponse);
			}
		}
	
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
			while(this.delegate.hasNext()) {
				final SAMRecord rec = delegate.next();
				final int rec_tid= rec.getReferenceIndex();
				if(rec_tid < this.interval.referenceIndex) {
					continue;
					}
				else if(rec_tid > this.interval.referenceIndex) {
					close();
					return null;
					}
				else if(rec.getStart() > this.interval.start) {
					close();
					return null;
					}
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
				LOG.debug("getting header for "+url);

				try( CloseableHttpResponse httpResponse = httpClient.execute(httpGet)) {
					 final int responseCode = httpResponse.getStatusLine().getStatusCode();
					LOG.debug("responseCode: "+responseCode);

					 if(responseCode != 200)
					 	{
						throw new RuntimeIOException("Response code was not 200. Detected response was "+responseCode);
					 	}
					
					 final HttpEntity entity = httpResponse.getEntity();
					 LOG.info("content length 0");
					 this.content_length = entity.getContentLength();
					 LOG.info("content length 1 "+this.content_length);
					
					 
					SamReader sr = null;
					    {
						sr= srf.open(SamInputResource.of(entity.getContent())); 
						LOG.info("get sr");
						this.samFileHeader = sr.getFileHeader();
						LOG.info("got sr");
						//in.close();
						//sr.close();
						}
					
					LOG.info("getFileHeader"+this.samFileHeader.getSortOrder());
				 
				// end try with resource
				if(!this.samFileHeader.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
					throw new IOException("Bam "+url+" is not sorted on coordinate: "+this.samFileHeader.getSortOrder());
					}
				}
			}
		
		public SAMFileHeader getFileHeader() {
			return samFileHeader;
		}

		
		private CloseableIterator<SAMRecord> queryAtOffset(final long start_offset) throws IOException {
			LOG.warn("opening at offset="+start_offset);

			if(start_offset>=this.content_length) return EOF_ITER;
			LOG.warn("read at "+start_offset);
			final byte buffer[]  =  new byte[byte_buffer_size];
            httpGet.setHeader("Range", "bytes="+start_offset+"-"+this.content_length);
            CloseableHttpResponse httpResponse = httpClient.execute(httpGet);
            final int responseCode = httpResponse.getStatusLine().getStatusCode();

            if(!(responseCode == 200 || responseCode == 206 /* partial */))
 			 	{
            	httpResponse.close();
            	throw new RuntimeIOException("Response code was not 200/206 for offset '"+start_offset+"'. Detected response was "+responseCode);
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
			 for(int offset=0;offset +1 < buffer.length;offset++) {
		         final ByteBuffer byteBuffer = ByteBuffer.wrap(buffer,offset,2);
		         byteBuffer.order(ByteOrder.LITTLE_ENDIAN);
		         if (byteBuffer.get() != BlockCompressedStreamConstants.GZIP_ID1) continue;
		         if (byteBuffer.get() != (byte)BlockCompressedStreamConstants.GZIP_ID2) continue;
		         
				 // https://twitter.com/jomarnz/status/1205533685413486594
				 mySamRecordIter = null;
				 final ByteArrayInputStream bais = new ByteArrayInputStream(buffer, offset, buffer.length);
				 final SequenceInputStream mergeIn= new SequenceInputStream(bais,in);
				 final BlockCompressedInputStream bci = new BlockCompressedInputStream(mergeIn);
				 try {
					 LOG.info("ok1");
					mySamRecordIter = new MySamRecordIterator(this.samFileHeader,bci,this.url);
					 LOG.info("ok2");
					mySamRecordIter.hasNext();//force checking next
					mySamRecordIter.onClose=()->{
						try {
							//in.close();
							httpResponse.close();
							}
						catch(final IOException err) 
							{
							LOG.warning(err);
							}
						};
					LOG.info("got it !");
					return mySamRecordIter;
				 	}
				 catch(java.lang.OutOfMemoryError err) {
					 mySamRecordIter = null;
					 err.printStackTrace();
				 	}
				 catch(final Throwable err) {
					 System.err.println("offset:"+offset+" "+err.getMessage()+" "+err.getClass());
					 mySamRecordIter = null;
				 	}
			 	}
			//in.close();
			httpResponse.close();
			LOG.info("bgzf signature no found");
			return EOF_ITER;
			}
		
		
		
		public CloseableIterator<SAMRecord> query(final Locatable locatable) throws IOException {
			LOG.info("ici "+locatable);
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(samFileHeader);
			final int tid = dict.getSequenceIndex(locatable.getContig());
			if(tid==-1) {
				LOG.info("unknown contig "+locatable.getContig());
				return EOF_ITER;
			}
			final QueryInterval userQueryInterval = new QueryInterval(tid, locatable.getStart(),locatable.getEnd());
			long byte_start = 0L;
			long byte_end=this.content_length;
			int repeat_dichotomy = dichotomy_repeat;
		    long len = byte_end - byte_start;
		    while (len > 0L)
		            {
		    		LOG.info("len="+len);
		    		repeat_dichotomy--;
		            final long half = len / 2;
		            final long middle = byte_start + half;
		            if(middle < 1_000_000L) {
		            	LOG.info("return default iter");
		            	return new MySamFilterIterator(new HeadingSamRecordIterator(this.httpGet),userQueryInterval);
		            }
		            
		            
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

		if(do_debug || true) {
			// https://hc.apache.org/httpcomponents-client-4.5.x/logging.html
			System.setProperty("org.apache.commons.logging.Log","org.apache.commons.logging.impl.SimpleLog");
			System.setProperty("org.apache.commons.logging.simplelog.showdatetime","true");
			System.setProperty("org.apache.commons.logging.simplelog.log.org.apache.http","DEBUG");
			System.setProperty("org.apache.commons.logging.simplelog.log.org.apache.http.wire","ERROR");
		}
		

		HttpHost proxy = new HttpHost("cache.ha.univ-nantes.fr", 3128, "http");
		/** create http client */
		this.httpClient = HttpClientBuilder.
				create().
				setUserAgent(IOUtils.getDefaultUserAgent()).
				build();
		
		
		try( CustomSamReader sr =  new CustomSamReader(urlStr)) {
			LOG.info("A");
			final SAMFileHeader samFileHeader = sr.getFileHeader();
			LOG.info("B");
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(samFileHeader);
			final Locatable userInterval = IntervalParserFactory.newInstance(dict).make().apply(this.intervalStr).orElse(null);
			if(userInterval==null) {
				LOG.error("cannot parse interval "+this.intervalStr+" for "+urlStr);
				return -1;
				}
			LOG.info("query "+userInterval);
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
