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
import java.net.URL;
import java.net.URLDecoder;
import java.util.Arrays;
import java.util.stream.Collectors;

import org.apache.http.Header;
import org.apache.http.HttpStatus;
import org.apache.http.auth.AuthScope;
import org.apache.http.auth.UsernamePasswordCredentials;
import org.apache.http.client.config.RequestConfig;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.methods.HttpHead;
import org.apache.http.client.methods.HttpRequestBase;
import org.apache.http.impl.client.BasicCredentialsProvider;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClientBuilder;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.client.protocol.HttpClientContext;

import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.seekablestream.ISeekableStreamFactory;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.CloserUtil;

/**
 * 
 * Custom ISeekableStreamFactory, handle user/password for URLs.
 * Used for BBFile constructor for remote data or HicFileReader,
 *
 */
public class CustomSeekableStreamFactory 
	implements ISeekableStreamFactory {
	private final ISeekableStreamFactory defaultInstance;
	private String user = null;
	private String password = null;
	// https://stackoverflow.com/questions/56598349
	private boolean normalizeURI = true;
	// user agent
	private String userAgent = null;
	// use httpGet instead of httpHead to get the content-length
	private boolean usingHttpHead = true;
	
	
	
	/** construct with htsjdk.samtools.seekablestream.SeekableStreamFactory.getInstance()</code> */
	public CustomSeekableStreamFactory()
		{
		this.defaultInstance = htsjdk.samtools.seekablestream.SeekableStreamFactory.getInstance();
		}
	
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
	
	/** https://stackoverflow.com/questions/56598349 default is TRUE */
	public CustomSeekableStreamFactory setNormalizeURI(boolean normalizeURI) {
		this.normalizeURI = normalizeURI;
		return this;
		}
	
	public boolean isNormalizeURI() {
		return normalizeURI;
		}
	
	/** amazon doesn't like httpHead https://stackoverflow.com/questions/56598349/wget-returns-200-ok-while-org-apache-http-impl-client-closeablehttpclient-return#comment99824950_56605768 */
	public CustomSeekableStreamFactory setUsingHttpHead(boolean usingHttpHead) {
		this.usingHttpHead = usingHttpHead;
		return this;
	}
	
	public boolean isUsingHttpHead() {
		return usingHttpHead;
	}
	
	public CustomSeekableStreamFactory setUserAgent(final String userAgent) {
		this.userAgent = userAgent;
		return this;
		}
	
	public String getUserAgent() {
		return userAgent;
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
		
		hb.useSystemProperties();
		
		if(!StringUtils.isBlank(getUserAgent())) {
			hb.setUserAgent(getUserAgent());
		}
		
		// https://stackoverflow.com/questions/56598349
		final RequestConfig requestConfig = RequestConfig.
				custom().
				setNormalizeUri(this.isNormalizeURI()).
				build();
		final HttpClientContext clientContext = HttpClientContext.create();
		clientContext.setRequestConfig(requestConfig);
		
		// set login/password
		if(!StringUtils.isBlank(p_user) && !StringUtils.isBlank(p_password)) {
			final BasicCredentialsProvider provider = new BasicCredentialsProvider();
			provider.setCredentials( AuthScope.ANY, new UsernamePasswordCredentials(p_user,p_password));
			hb.setDefaultCredentialsProvider(provider);
			}
			
		
		final CloseableHttpClient httpClient = hb.build();
		
		final HttpRequestBase httpHead = this.isUsingHttpHead()?
				  new HttpHead(url.toExternalForm())
				: new HttpGet(url.toExternalForm())
				;
		try
			{
			final CloseableHttpResponse response = httpClient.execute(httpHead,clientContext);
			if (response.getStatusLine().getStatusCode() != HttpStatus.SC_OK ) {
				if(response.getAllHeaders()!=null)
					{
					/*
					System.err.println( Arrays.stream(response.getAllHeaders()).
							map(X->X.getName()+":"+X.getValue()).
							collect(Collectors.joining(" ; ")));
					*/
					}
				
				final String msg = "Unexpected Http status code "
     		            + response.getStatusLine()+" for \""+ url+"\""
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
			
			final Header acceptRanges[]  = response.getHeaders("Accept-Ranges");
			if(acceptRanges==null || Arrays.stream(acceptRanges).noneMatch(S->S.getValue().equals("bytes"))) {
				throw new IOException("Cannot get Accept-Ranges: bytes "+
						acceptRanges==null?"": Arrays.stream(acceptRanges).
								map(X->X.getName()+":"+X.getValue()).
								collect(Collectors.joining(" ; "))
						);
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
				httpClient,
				clientContext
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
