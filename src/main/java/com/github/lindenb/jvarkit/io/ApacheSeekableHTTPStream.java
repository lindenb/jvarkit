/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.io;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

import org.apache.http.HttpEntity;
import org.apache.http.HttpStatus;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.client.protocol.HttpClientContext;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.CloserUtil;


/**
 * an implentation of SeekableStream using apache http-components
 *
 */
class ApacheSeekableHTTPStream extends SeekableStream {
	private static final int BUFFER_SIZE= 2_000_000;
	private final byte buffer[]=new byte [BUFFER_SIZE];
	private long buffer_offset = 0;
	private int buffer_size = 0;
	private int buffer_index = BUFFER_SIZE;
	
    // private long position = 0L;
    private final long contentLength ;
    private final URL url;
    private final CloseableHttpClient httpClient;
    // https://stackoverflow.com/questions/56598349/
    private final HttpClientContext clientContext;

    ApacheSeekableHTTPStream(
    	final URL url,
    	final long contentLength,
    	final CloseableHttpClient httpClient,
    	final HttpClientContext clientContextOrNull
    	) {
        this.url = url;
        this.contentLength = contentLength;
        this.httpClient = httpClient;
        this.clientContext = (clientContextOrNull == null? HttpClientContext.create() : clientContextOrNull);
    	}

    @Override
    public long position() {
        return this.buffer_offset + this.buffer_index;
    }

    @Override
    public long length() {
        return this.contentLength;
    }

    @Override
    public long skip(final long n) throws IOException {
    	final long bytesToSkip = Math.min(n, this.contentLength - position());
        final long new_offset = this.buffer_offset + bytesToSkip;
       
        
        if(new_offset >= this.buffer_offset && new_offset < this.buffer_offset + this.buffer_size) {
    		this.buffer_index = (int)(new_offset-this.buffer_offset);
    		}
    	else {
	        this.buffer_offset = new_offset;
	        this.buffer_index = 0;
	        this.buffer_size = 0;
	    	}
        
        return bytesToSkip;
    }

    @Override
    public boolean eof() throws IOException {
        return this.contentLength > 0 && this.position() >= this.contentLength;
    }

    @Override
    public void seek(final long position) {
    	if(position >= this.buffer_offset && position < this.buffer_offset + this.buffer_size) {
    		this.buffer_index = (int)(position-this.buffer_offset);
    		}
    	else {
	        this.buffer_offset = position;
	        this.buffer_index = 0;
	        this.buffer_size = 0;
	    	}
    	}
    
    
    
    private int readNextByte() throws IOException {
    	if(this.buffer_index < this.buffer_size) {
    		final int c = this.buffer[this.buffer_index];
    		this.buffer_index++;
    		return c  & 0xFF;
    		}
    	//System.err.println("refill from "+this.buffer_offset);
    	this.buffer_offset += this.buffer_size;
    	this.buffer_index = 0;
    	this.buffer_size = 0;
    	
    	final int max_n_bytes_to_reads = (int)Math.min((long)this.buffer.length,this.contentLength - this.buffer_offset);
    	
        if (max_n_bytes_to_reads == 0 ) return -1;
        if (this.buffer_offset == this.contentLength) return -1;

        CloseableHttpResponse httpResponse = null;
        InputStream is = null;
        
        try {
        	final HttpGet httpGet = new HttpGet(this.url.toExternalForm());
        	
        	

            long endRange = this.buffer_offset + max_n_bytes_to_reads - 1;
            // IF we know the total content length, limit the end range to that.
            if (contentLength > 0) {
                endRange = Math.min(endRange, contentLength);
            }
            final String byteRange = "bytes=" + this.buffer_offset + "-" + endRange;
            
            httpGet.addHeader("Range", byteRange);
          
            httpResponse = this.httpClient.execute(httpGet,this.clientContext);
            if (httpResponse.getStatusLine().getStatusCode() != HttpStatus.SC_PARTIAL_CONTENT ) {
				final String msg = "Unexpected Http status code "
     		            + httpResponse.getStatusLine()+" for "+ url +" in range "+byteRange;
				httpResponse.close();
				httpResponse = null;
 		        throw new IOException(msg);
 		   		}
            final HttpEntity entity = httpResponse.getEntity();
            
            is = entity.getContent();

            while (this.buffer_size < max_n_bytes_to_reads) {
                int count = is.read(this.buffer, this.buffer_size, max_n_bytes_to_reads - this.buffer_size);
                if (count < 0) {
                    if (this.buffer_size == 0) {
                        return -1;
                    } else {
                        break;
                    }
                }
                this.buffer_size += count;
            	}

            int c= this.buffer[this.buffer_index];
            this.buffer_index++;
            
        	//System.err.println("refilled from "+this.buffer_offset+" size="+this.buffer_size);

            return c  & 0xFF;
        	}
        finally {
            CloserUtil.close(is);
            CloserUtil.close(httpResponse);
        }
    
    	
    }
    

    @Override
    public int read(byte[] array, int offset, int len) throws IOException {
    	int n = 0;
        while (n < len) {
        	int c = this.readNextByte();
        	if(c==-1) {
    		  if (n == 0) {
                  return -1;
              } else {
                  break;
              }
        	}
        	
        	array[ offset+n ] = (byte)c;
            n++;
        	}
        return n;
    	}
   

    @Override
    public void close() throws IOException {
        CloserUtil.close(this.httpClient);
    }


    @Override
    public int read() throws IOException {
    	return readNextByte();
    }

    @Override
    public String getSource() {
        return this.url.toString();
    }
    
    
   }
