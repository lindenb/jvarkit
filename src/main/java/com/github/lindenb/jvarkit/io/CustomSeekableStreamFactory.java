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
import java.net.URL;
import java.net.URLDecoder;

import org.apache.http.Header;
import org.apache.http.HttpStatus;
import org.apache.http.auth.AuthScope;
import org.apache.http.auth.UsernamePasswordCredentials;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpHead;
import org.apache.http.impl.client.BasicCredentialsProvider;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClientBuilder;
import org.apache.http.impl.client.HttpClients;

import htsjdk.samtools.seekablestream.ISeekableStreamFactory;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;

/**
 * 
 * Custom ISeekableStreamFactory, handle user/password for URLs.
 * Used for BBFile constructor for remote data
 *
 */
public class CustomSeekableStreamFactory 
	implements ISeekableStreamFactory {
	private final ISeekableStreamFactory defaultInstance;
	private String user = null;
	private String password = null;

	/** construct with something like <code>new CustomSeekableStreamFactory(htsjdk.samtools.seekablestream.SeekableStreamFactory.getInstance());</code> */
	public CustomSeekableStreamFactory(final ISeekableStreamFactory defaultInstance)
		{
		this.defaultInstance = defaultInstance;
		if(defaultInstance==null) throw new IllegalArgumentException("default instance cannot be null");
		}
	/** get user for  http authentication */
	public String getUser() {
		return user;
		}
	/** get password for http authentication */
	public String getPassword() {
		return password;
		}
	
	/** set user for  http authentication */
	public CustomSeekableStreamFactory setUser(final String u) {
		this.user = u;
		return this;
		}
	
	/** set password for http authentication */
	public CustomSeekableStreamFactory setPassword(final String p) {
		this.password = p;
		return this;
		}

	/** reset user/password */
	public CustomSeekableStreamFactory reset() {
		return setUser(null).setPassword(null);
		}
	
	private ISeekableStreamFactory getDelegate() {
		return this.defaultInstance;
	}
	
    @Override
    public SeekableStream getStreamFor(final URL url) throws IOException {
    	return this.getStreamFor(url.toExternalForm());//don't call delegate
    	}

    @Override
    public SeekableStream getStreamFor(final String path) throws IOException {
    	if (path.startsWith("http:") || path.startsWith("https:")) {
    		 String p_user=null,p_password=null;
    		 final URL url;
    		 int slashslash = path.indexOf("://");
    		 int colon =  slashslash==-1?-1:path.indexOf(':',slashslash+3);
    		 int arobase = colon==-1?-1:path.indexOf('@',colon+1);
    		 int dot = colon==-1?-1:path.indexOf('.',arobase+1);
    		 int firstslash =  slashslash==-1?-1:path.indexOf('/',slashslash+3);
    		 if(arobase!=-1 && colon!=-1 && slashslash<colon && colon< arobase && arobase<dot && dot < firstslash)
    		 	{
    			url=new URL(path.substring(0, slashslash+3)+path.substring(arobase+1));
    			p_user = URLDecoder.decode(path.substring(slashslash+3, colon), "UTF-8");
    			p_password = URLDecoder.decode(path.substring(colon+1, arobase), "UTF-8");
    		 	}
    		 else
    		 	{
    			url=new URL(path);
    			p_user = this.getUser();
    			p_password = this.getPassword();
    		 	}
    		if(StringUtil.isBlank(p_user) || StringUtil.isBlank(p_password)) 
    			{
    			return  getDelegate().getStreamFor(path);
    			}
    		 return openHttp(url,p_user,p_password);
	    	 }
    	return getDelegate().getStreamFor(path);
    	}
    
    @Override
    public SeekableStream getBufferedStream(final SeekableStream stream){
    	return getDelegate().getBufferedStream(stream);
    }

    @Override
    public SeekableStream getBufferedStream(final SeekableStream stream, int bufferSize){
    	return getDelegate().getBufferedStream(stream,bufferSize);
    }

    private SeekableStream openHttp(final URL url,final String p_user,final String p_password) 
    	throws IOException {
		final HttpClientBuilder hb = HttpClients.custom();
		final BasicCredentialsProvider provider = new BasicCredentialsProvider();
		provider.setCredentials( AuthScope.ANY, new UsernamePasswordCredentials(p_user,p_password));
		hb.setDefaultCredentialsProvider(provider);
			
		
		final CloseableHttpClient httpClient = hb.build();
		
		final HttpHead httpHead = new HttpHead(url.toExternalForm());
		try
			{
			final CloseableHttpResponse response = httpClient.execute(httpHead);
			if (response.getStatusLine().getStatusCode() != HttpStatus.SC_OK ) {
				final String msg = "Unexpected Http status code "
     		            + response.getStatusLine()+" for "+ url
     		            ;
				response.close();
				httpClient.close();
 		        throw new IOException(msg);
 		   		}
			
			final Header contentLengthHeader = response.getFirstHeader("Content-Length");
			if(contentLengthHeader==null)
				{
				final String msg = "Cannot get Content-length for "+ url;
				response.close();
				httpClient.close();
 		        throw new IOException(msg);
				}
			
			long contentLength;
			try {
				contentLength = Long.parseLong(contentLengthHeader.getValue());
				if(contentLength<0) throw new NumberFormatException("Negative content length for "+url);
				}
			catch(final NumberFormatException err)
				{
				final String msg = "bad Content-Length in "+contentLengthHeader+" for "+ url;
				response.close();
				httpClient.close();
 		        throw new IOException(msg,err);
				}
			response.close();
			final ApacheSeekableHTTPStream stream = new ApacheSeekableHTTPStream(
				url,
				contentLength,
				httpClient
				);
			return stream;
			}
		catch(final IOException err)
			{
			CloserUtil.close(httpClient);
			throw err;
			}		
    	}
    
}
