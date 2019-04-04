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
package com.github.lindenb.jvarkit.util;

import java.io.InputStream;
import java.net.URL;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Set;
import java.util.jar.Attributes;
import java.util.jar.Manifest;

import com.github.lindenb.jvarkit.util.jcommander.Launcher;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

/** get general data about jvarkit in java manifest */
public class JVarkitVersion {
private static JVarkitVersion INSTANCE=null;
private String compilationDate = null;
private String gitHash = null;
private String htsjdkVersion = null;

private JVarkitVersion() {
	Enumeration<URL> resources = null;
	try {
		resources = getClass().getClassLoader().getResources("META-INF/MANIFEST.MF");
		}
	catch(Exception err) {
		resources = null;
		}
	while (resources!=null && resources.hasMoreElements()) {
		InputStream is = null;
	    try {
	        is = resources.nextElement().openStream();
	        final Manifest manifest = new Manifest(is);
	        final  Attributes attr = manifest.getMainAttributes();
	        if(StringUtil.isBlank(this.compilationDate))
	        	{
	        	final String cp=attr.getValue("Compile-Date");
		        if(!StringUtil.isBlank(cp))
		        	{
		        	compilationDate = cp;
		        	}
	        	}
	        if(StringUtil.isBlank(this.gitHash))
	        	{
	        	final String cp=attr.getValue("Git-Hash");
		        if(!StringUtil.isBlank(cp))
		        	{
		        	gitHash = cp;
		        	}
	        	}
	        if(StringUtil.isBlank(this.htsjdkVersion))
	        	{
	        	final String cp=attr.getValue("Htsjdk-Version");
		        if(!StringUtil.isBlank(cp))
		        	{
		        	htsjdkVersion = cp;
		        	}
	        	}
	    	} 
	    catch (final Exception err) {
	    	//
	    	}
	    finally
	    	{
	    	CloserUtil.close(is);
	    	}
		}
	if(StringUtil.isBlank(this.compilationDate)) {
		compilationDate = "(undefined)";
		}
	if(StringUtil.isBlank(this.gitHash)) {
		gitHash = "(undefined)";
		}
	if(StringUtil.isBlank(this.htsjdkVersion)) {
		htsjdkVersion = "(undefined)";
		}
	}

/** return compilation date, never null */
public String getCompilationDate() {
	return compilationDate;
	}

/** return git hash , never null */
public String getGitHash() {
	return gitHash;
	}
/** return htsjdk version , never null */
public String getHtsjdkVersion() {
	return htsjdkVersion;
	}

public static JVarkitVersion getInstance() {
	if(INSTANCE==null) {
		synchronized (JVarkitVersion.class) {
			if(INSTANCE==null) {
				INSTANCE = new JVarkitVersion();
				}
			}
		}
	return INSTANCE;
	}

public Set<VCFHeaderLine> getMetaData(String prefix) {
	if(StringUtil.isBlank(prefix)) prefix="jvarkit";
	if(!prefix.endsWith(".")) prefix+=".";
	final Set<VCFHeaderLine> metaData = new HashSet<>();
	metaData.add(new VCFHeaderLine(prefix+"meta",getLabel()));
	return metaData;
	}

public VCFHeader addMetaData(String prefix,final VCFHeader header) {
	getMetaData(prefix).stream().forEach(H->header.addMetaDataLine(H));
	return header;
	}
public VCFHeader addMetaData(final Launcher app,final VCFHeader header) {
	header.addMetaDataLine(new VCFHeaderLine(
			app.getProgramName()+".meta",
			getLabel()+" cmd:"+
			app.getProgramCommandLine().replaceAll("[\\\n\r\"\']+", " ").trim())
			);
	return header;
	}

public SAMFileHeader addMetaData(final Launcher app,final SAMFileHeader header) {
	header.addComment(
			app.getProgramName()+". "+
			getLabel()+". cmd:"+
			app.getProgramCommandLine().replaceAll("[\\\n\r\"\']+", " ").trim()
			);
	return header;
	}


public String getLabel() {
	return String.join(" ",
			"compilation:"+ getCompilationDate(),
			"githash:"+ getGitHash(),
			"htsjdk:"+ getHtsjdkVersion(),
			"date:" + new SimpleDateFormat("yyyyMMddHHmmss").format(new Date())
			);
	
	
}

@Override
public String toString() { return getLabel();}
}
