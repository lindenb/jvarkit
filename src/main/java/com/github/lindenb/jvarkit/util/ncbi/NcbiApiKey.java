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
package com.github.lindenb.jvarkit.util.ncbi;

import java.io.File;
import java.io.FileReader;
import java.util.Properties;
import java.util.regex.Pattern;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
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
	public static final String ENV_NAME="NCBI_API_KEY";
	public static final String JVM_PARAM="ncbi.api.key";
	private static final Logger LOG=Logger.build(NcbiApiKey.class).make();
	@Parameter(names="--ncbi-api-key",description=
			"NCBI API Key see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ ."
			+ 	"If undefined, it will try to get in that order: "
			+ " 1) environment variable ${"+ENV_NAME+"} ; "
			+ " 2) the jvm property \""+JVM_PARAM+"\" ;"
			+ "	3) A java property file ${HOME}/"+CONFIG_FILE+" and key "+PARAM
			)
	private String key = null;
	private boolean _searched = false;
	
	public NcbiApiKey()
		{
		
		}
	
	public NcbiApiKey(final String s)
		{
		this.key = s;
		}
	
	
	public String getApiKey() {		
		if(!StringUtil.isBlank(this.key)) return this.key;
		if(this._searched) return null;
		this._searched = true;
		this.key = System.getenv(ENV_NAME);
		if(!StringUtil.isBlank(this.key))
			{
			return this.key;
			}
		this.key = System.getProperty(JVM_PARAM,null);
		if(!StringUtil.isBlank(this.key))
			{
			return this.key;
			}
		final File keyFile = new File(System.getProperty("user.home"),CONFIG_FILE);
		FileReader r=null;
		try {
			if(keyFile.exists()) {
				final Properties prop = new Properties();
				r=new FileReader(keyFile);
				prop.load(r);
				r.close();r=null;
				this.key=prop.getProperty(PARAM,null);
				if(!StringUtil.isBlank(this.key))
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
		LOG.warn("\n"+
			"*****\n"+
			"*\n"+
			"* NCBI api_key is undefined. see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/.\n"+
			"* 1) can find it in environment variable  `export "+ENV_NAME+"=xxxxxxxxx`\n" +
			"* 2) can find it in jvm property java -D"+JVM_PARAM+"=xxxxxxxxx ....\n" +
			"* 3) can find it in java property file ( https://en.wikipedia.org/wiki/.properties ) "+keyFile+" with property "+PARAM+"=xxxxxxxxx\n" +
			"*\n"+
			"*****"
			);
		
		
		return this.key;	
		}
	
	/** throws an exeception if api key is not defined*/
	public void assertApiKeyDefined()
		{
		if(!isApiKeyDefined())
			{
			throw new JvarkitException.UserError(
					"NCBI API Key is undefined. see parameters and https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ ");
			}
		}
	
	/** return true if api key is defined and not empty/null */
	public boolean isApiKeyDefined()
		{
		return !StringUtil.isBlank(this.getApiKey());
		}
	
	/** return '&api_key=xxxxx' or empty string if no key is defined */
	public String getAmpParamValue() {
		final String s=this.getApiKey();
		if(StringUtil.isBlank(s)) return "";
		return "&"+PARAM+"="+s;
	}
	
	/** hide the api key if it is present in a string */
	public String mask(final String s) {
		if(s==null ) return s;
		final String kv = getApiKey();
		if(StringUtil.isBlank(kv)) return s;
		return s.replaceAll(Pattern.quote(kv),"xxxxxxx");
		}
}
