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
package com.github.lindenb.jvarkit.io;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

import org.apache.http.HttpEntity;
import org.apache.http.HttpStatus;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.CloserUtil;

public class ApacheSeekableHTTPStream extends SeekableStream {

    private long position = 0L;
    private final long contentLength ;
    private final URL url;
    private final CloseableHttpClient httpClient;

    /* use CustomSeekableStreamFactory */ ApacheSeekableHTTPStream(
    	final URL url,
    	final long contentLength,
    	final CloseableHttpClient httpClient
    	) {
        this.url = url;
        this.contentLength = contentLength;
        this.httpClient = httpClient;
    }

    @Override
    public long position() {
        return position;
    }

    @Override
    public long length() {
        return contentLength;
    }

    @Override
    public long skip(final long n) throws IOException {
    	final long bytesToSkip = Math.min(n, contentLength - position);
        this.position += bytesToSkip;
        return bytesToSkip;
    }

    @Override
    public boolean eof() throws IOException {
        return this.contentLength > 0 && this.position >= this.contentLength;
    }

    @Override
    public void seek(final long position) {
        this.position = position;
    }

    @Override
    public int read(byte[] buffer, int offset, int len) throws IOException {

        if (offset < 0 || len < 0 || (offset + len) > buffer.length) {
            throw new IndexOutOfBoundsException("Offset="+offset+",len="+len+",buflen="+buffer.length);
        }
        if (len == 0 ) {
            return 0;
        }
        if (this.position == this.contentLength) {
            return -1;
        }

        CloseableHttpResponse httpResponse = null;
        InputStream is = null;
        int n = 0;
        try {
        	final HttpGet httpGet = new HttpGet(this.url.toExternalForm());
        	
        	

            long endRange = position + len - 1;
            // IF we know the total content length, limit the end range to that.
            if (contentLength > 0) {
                endRange = Math.min(endRange, contentLength);
            }
            final String byteRange = "bytes=" + position + "-" + endRange;
            
            httpGet.addHeader("Range", byteRange);
          
            httpResponse = this.httpClient.execute(httpGet);
            if (httpResponse.getStatusLine().getStatusCode() != HttpStatus.SC_PARTIAL_CONTENT ) {
				final String msg = "Unexpected Http status code "
     		            + httpResponse.getStatusLine()+" for "+ url +" in range "+byteRange;
				httpResponse.close();
				httpResponse = null;
 		        throw new IOException(msg);
 		   		}
            final HttpEntity entity = httpResponse.getEntity();
            
            is = entity.getContent();

            while (n < len) {
                int count = is.read(buffer, offset + n, len - n);
                if (count < 0) {
                    if (n == 0) {
                        return -1;
                    } else {
                        break;
                    }
                }
                n += count;
            	}
            this.position += n;

            return n;
        	}
        finally {
            CloserUtil.close(is);
            CloserUtil.close(httpResponse);
        }
    }


    @Override
    public void close() throws IOException {
        CloserUtil.close(this.httpClient);
    }


    @Override
    public int read() throws IOException {
    	byte []tmp=new byte[1];
    	read(tmp,0,1);
    	return (int) tmp[0] & 0xFF; 
    }

    @Override
    public String getSource() {
        return this.url.toString();
    }
    
    
   }
