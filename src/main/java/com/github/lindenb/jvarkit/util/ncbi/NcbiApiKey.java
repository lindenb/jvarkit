/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.util.ncbi;

import java.io.File;
import java.io.FileReader;
import java.util.Properties;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;

/**
see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/?campaign=twitter-11022017

> If you regularly use the E-utilities API, we have important news 
> for you: NCBI is now providing API keys for the E-utilities! 
> After May 1, 2018, NCBI will limit your access to the E-utilities
> unless you have one of these keys. Obtaining an API key is quick,
> and simple, and will allow you to access NCBI data faster.
> If you don’t have an API key, E-utilities will still work, but
> you may be limited to fewer requests than allowed with an API key.


 */
public class NcbiApiKey {
	private static final String CONFIG_FILE=".ncbi.properties";
	public static final String PARAM="api_key";
	private static final Logger LOG=Logger.build(NcbiApiKey.class).make();
	@Parameter(names="--ncbi-api-key",description="NCBI API Key see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ . If undefined. Will try to read in that order: 1) A java XML property file ${HOME}/"+CONFIG_FILE+" and key "+PARAM+" 2) the jvm property \"ncbi.api.key\" 3) environment variable NCBI_API_KEY")
	private String key = null;
	private boolean _searched = false;
	
	public String getApiKey() {		
		if(this.key!=null) return this.key;
		if(this._searched) return null;
		this._searched = true;
		FileReader r=null;
		try {
			final File keyFile = new File(System.getProperty("user.home, def"),CONFIG_FILE);
			if(keyFile.exists()) {
				final Properties prop = new Properties();
				r=new FileReader(keyFile);
				prop.load(r);
				r.close();r=null;
				if(!prop.containsKey(PARAM)) {
					LOG.warning("property \""+PARAM+"\" not defined in "+keyFile);
					}
				this.key=prop.getProperty(PARAM);
				if(StringUtil.isBlank(this.key))
					{
					this.key=null;
					LOG.warning("property \""+PARAM+"\" is empty in "+keyFile);
					}
				else
					{
					return this.key;
					}
				}
			else
				{
				LOG.warn("NCBI property file \""+keyFile+"\" not found.");
				}	
			}
		catch(final Exception err) {
			LOG.error(err);
			}
		finally
			{
			CloserUtil.close(r);
			}
		this.key = System.getProperty("ncbi.api.key");
		if(StringUtil.isBlank(this.key))
			{
			LOG.warning("JVM property \"ncbi.api.key\" is empty");
			}
		else
			{
			return this.key;
			}
		this.key = System.getenv("NCBI_API_KEY");
		if(StringUtil.isBlank(this.key))
			{
			LOG.warning("env property \"NCBI_API_KEY\" is empty");
			}
		return this.key;	
		}
	/** return '&api_key=xxxxx' or empty string if no key is defined */
	public String getAmpParamValue() {
		final String s=this.getApiKey();
		if(StringUtil.isBlank(s)) return "";
		return "&"+PARAM+"="+s;
	}
}
