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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.ngsfiles;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.xml.stream.XMLEventFactory;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamResult;

import htsjdk.variant.vcf.VCFHeader;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.illumina.FastQName;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;
import com.sleepycat.bind.tuple.StringBinding;
import com.sleepycat.je.Cursor;
import com.sleepycat.je.Database;
import com.sleepycat.je.DatabaseConfig;
import com.sleepycat.je.DatabaseEntry;
import com.sleepycat.je.Environment;
import com.sleepycat.je.EnvironmentConfig;
import com.sleepycat.je.LockMode;
import com.sleepycat.je.OperationStatus;
import com.sleepycat.je.Transaction;
/*

BEGIN_DOC

## Example

### Example 01 : scanning all files:

we scan all files under /common/data/projects.
all missing files/dir previously inserted will be deleted from the berkeley database

```bash
$ java -jar dist/ngsfilesscanner.jar -B /path/to/bdbdir /commun/data/projects/
```
### Example 02 : dumping the database

we scan all files under /common/data/projects.

```bash
$ java -jar dist/ngsfilesscanner.jar -B /path/to/bdbdir -D 
```
```xml
<?xml version="1.0" encoding="UTF-8"?>
<ngs-files>
  <fastq-dir directory="/commun/data/projects/Project_P1/Samplex1">
    <samples>
      <sample size="2834419615">x1</sample>
    </samples>
  </fastq-dir>
  (...)
 <vcf timestamp="1398643093000" file="/commun/data/projects/path/Samples/S2/S2.varscan.annotations.vcf.gz" filename="S2.varscan.annotations.vcf.gz" modified="Mon Apr 28 01:58:13 CEST 2014" size="21053412">
    <samples>
      <sample>S2</sample>
    </samples>
  </vcf>
</ngs-files>
```
## setting-up a CRON job

content of `/etc/cron.daily/ngsfilesscanner.cron`
```bash
#!/bin/bash
/usr/bin/java -cp /path/to/picard-tools-1.100/picard-1.100.jar:/path/to/picard-tools-1.100/sam-1.100.jar:/path/to/picard-tools-1.100/tribble-1.100.jar:/path/to/picard-tools-1.100/variant-1.100.jar:/path/to/picard-tools-1.100/commons-jexl-2.1.1.jar:/path/to/picard-tools-1.100/commons-logging-1.1.1.jar:/path/to/je-5.0.34/lib/je-5.0.34.jar:/path/to/jvarkit-git/ngsfilesscanner.jar com.github.lindenb.jvarkit.tools.ngsfiles.NgsFilesScanner -B /var/www/cgi-bin/ngsfiles /commun/data/
```
set permissions
```
$ sudo chown root:root /etc/cron.daily/ngsfilesscanner.cron
$ sudo chmod u+x ngsfilesscanner.cron
```
## creating a CGI dumping the results:

content of `/var/www/cgi-bin/ngsfiles.cgi`
```bash
#!/bin/bash

echo "Content-Type: text/xml"
echo ""
/usr/bin/java -cp /path/to/picard-tools-1.100/picard-1.100.jar:/path/to/picard-tools-1.100/sam-1.100.jar:/path/to/picard-tools-1.100/tribble-1.100.jar:/path/to/picard-tools-1.100/variant-1.100.jar:/path/to/picard-tools-1.100/commons-jexl-2.1.1.jar:/path/to/picard-tools-1.100/commons-logging-1.1.1.jar:/path/to/je-5.0.34/lib/je-5.0.34.jar:/path/to/jvarkit-git/ngsfilesscanner.jar com.github.lindenb.jvarkit.tools.ngsfiles.NgsFilesScanner -B /var/www/cgi-bin/ngsfiles  -D
```

set permissions
```
$ sudo chown root:root /var/www/cgi-bin/ngsfiles.cgi
$ sudo chmod u+x /var/www/cgi-bin/ngsfiles.cgi
```

check URL: `curl -s http://localhost/cgi-bin/ngsfiles.cgi | xmllint --format -`
```xml
<?xml version="1.0" encoding="UTF-8"?>
<ngs-files>
  <fastq-dir directory="/commun/data/projects/Project_P1/Samplex1">
    <samples>
      <sample size="2834419615">x1</sample>
    </samples>
  </fastq-dir>
  (...)
 <vcf timestamp="1398643093000" file="/commun/data/projects/path/Samples/S2/S2.varscan.annotations.vcf.gz" filename="S2.varscan.annotations.vcf.gz" modified="Mon Apr 28 01:58:13 CEST 2014" size="21053412">
    <samples>
      <sample>S2</sample>
    </samples>
  </vcf>
</ngs-files>
```

END_DOC
 
 */
@Program(name="ngsfilesscanner",
	description="Build a persistent database of NGS file. Dump as XML.",
	keywords={"ngs","bam","sam","vcf","xml"}
)
public class NgsFilesScanner extends AbstractScanNgsFilesProgram
	{
	private static final Logger LOG = Logger.build(NgsFilesScanner.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	
	private static final String SUFFIXES[]=new String[]{".bam",".vcf",".vcf.gz"};
	private Transaction txn=null;
	private Database database=null;
	private FileFilter fileFilter=new FileFilter()
		{
		@Override
		public boolean accept(File f) {
			if(f==null) return false;
			if(!f.canRead()) return false;
			if(f.isDirectory()) return true;
			String filename=f.getName().toLowerCase();
			if(filename.contains("jeter")) return false;
			if(filename.contains("_delete_")) return false;
			for(String suff: SUFFIXES)
				if(filename.endsWith(suff)) return true;
			return false;
			}
		};
		
		private FileFilter fastqFilter=new FileFilter()
			{
			@Override
			public boolean accept(File f) {
				if(f==null) return false;
				if(!f.canRead()) return false;
				if(!f.isFile()) return true;
				String filename=f.getName().toLowerCase();
				if(filename.contains("jeter")) return false;
				if(filename.contains("_delete_")) return false;
				if(filename.endsWith(".fastq.gz")) return true;
				if(filename.endsWith(".fastq")) return true;
				if(filename.endsWith(".fq.gz")) return true;
				if(filename.endsWith(".fqz")) return true;
				return false;
				}
			};	
		
	static final String DATABASE_NAME="ngsfile.db";
	
    private NgsFilesScanner()
    	{
    	
    	}		
    
    private void writeFile(XMLStreamWriter out,File f) throws XMLStreamException
    	{
    	out.writeAttribute("file", f.getAbsolutePath());
    	out.writeAttribute("filename", f.getName());
    	out.writeAttribute("size", String.valueOf(f.length()));
    	out.writeAttribute("modified", String.valueOf(new Date(f.lastModified())));
    	out.writeAttribute("timestamp", String.valueOf(f.lastModified()));
    	
    	}
    
    private void put(File f,String xml)
    	{
    	LOG.info("insert "+f);
		DatabaseEntry key=new DatabaseEntry();
		DatabaseEntry data=new DatabaseEntry();
		StringBinding.stringToEntry(f.getAbsolutePath(), key);
		StringBinding.stringToEntry(xml,data);
		this.database.put(this.txn, key, data);
    	}
    
    @Override
    protected void readFastq(File f) {
    	//do nothing
    	}
    
	@Override
    protected void readBam(final File f)
    	{
    	if(!f.canRead()) return;
    	SamReader r=null;
    	try {
    		StringWriter sw=new StringWriter();
    		XMLOutputFactory xof=XMLOutputFactory.newFactory();
    		XMLStreamWriter out=xof.createXMLStreamWriter(new StreamResult(sw));
    		
    		
    		out.writeStartElement("bam");
    		writeFile(out,f);
    		
			r=super.openSamReader(f.getPath());
			final SAMFileHeader h=r.getFileHeader();
			
			out.writeStartElement("samples");
			if(h!=null && h.getReadGroups()!=null)
				{
				Set<String> seen=new HashSet<String>();
				for(SAMReadGroupRecord rg: h.getReadGroups())
					{
					String sample=rg.getSample();
					if(sample==null || sample.isEmpty() || seen.contains(sample)) continue;
					
					seen.add(sample);
					out.writeStartElement("sample");
					out.writeCharacters(sample);
					out.writeEndElement();
					}
				}
			out.writeEndElement();
			
			out.writeEndElement();
			out.flush();
			out.close();
			sw.flush();
			
			put(f,sw.toString());
			} 
    	catch (Exception e)
    		{
    		LOG.warning(e);
			}
    	finally
    		{
    		CloserUtil.close(r);
    		}
    	}
   
   
    
   @Override
   protected void readVCF(File f)
		{
    	if(!f.canRead()) return;
    	LOG.debug("readVCF  "+f);
    	    	

    	VCFIterator r=null;
    	InputStream in=null;
    	try
    		{
    	
    		r=VCFUtils.createVCFIteratorFromFile(f);
        	VCFHeader header=r.getHeader();
        	
        	
    		StringWriter sw=new StringWriter();
    		XMLOutputFactory xof=XMLOutputFactory.newFactory();
    		XMLStreamWriter out=xof.createXMLStreamWriter(new StreamResult(sw));
    		out.writeStartElement("vcf");
    		writeFile(out,f);

    		out.writeStartElement("samples");
        	for(String sample:header.getSampleNamesInOrder())
	        	{
        		out.writeStartElement("sample");
				out.writeCharacters(sample);
				out.writeEndElement();
	    		}
        	out.writeEndElement();
        	
        	out.writeEndElement();
			out.flush();
			out.close();
			sw.flush();
			
			put(f,sw.toString());
    		}
    	catch(Exception err)
    		{
    		LOG.error(err);
    		}
    	finally
    		{
    		CloserUtil.close(r);
    		CloserUtil.close(in);
    		}
    	
		}
   
   
   private void fastqDir(File dir,Counter<String> samples)
   	{
   	try {
   		StringWriter sw=new StringWriter();
   		XMLOutputFactory xof=XMLOutputFactory.newFactory();
   		XMLStreamWriter out=xof.createXMLStreamWriter(new StreamResult(sw));
   		
   		
   		out.writeStartElement("fastq-dir");
   		out.writeAttribute("directory", dir.getAbsolutePath());
		out.writeStartElement("samples");
		for(String sample:samples.keySet())
			{
			out.writeStartElement("sample");
			out.writeAttribute("size",String.valueOf(samples.count(sample)));
			out.writeCharacters(sample);
			out.writeEndElement();
			}
			
		out.writeEndElement();
		
		out.writeEndElement();
		out.flush();
		out.close();
		sw.flush();
		
		put(dir,sw.toString());
		} 
   	catch (Exception e)
   		{
   		LOG.warning(e);
		}
   	finally
   		{
   		}
   	
   }
   
   private static final  Counter<String> EMPTY_COUNTER=new Counter<String>();
   
   private Counter<String> getFastqSampleInDirectory(File f)
	   	{
	   if(f==null || !f.isDirectory() || !f.exists() || !f.canRead()) return EMPTY_COUNTER;
	    File array[]=f.listFiles(this.fastqFilter);
	   	if(array==null) return EMPTY_COUNTER;
	   	Counter<String> fastqSamples=new Counter<String>();
	   	for(File f2:array)
	   		{
	   		FastQName fqName=FastQName.parse(f2);
			if(fqName.isValid() && fqName.getSample()!=null)
				{
				fastqSamples.incr(fqName.getSample(),f2.length());
				}
	   		}
	   	return fastqSamples;
	   	}
	
    private long lastDirTimeMillis=System.currentTimeMillis();
    private void recursive(File f) throws IOException
    	{
    	if(f==null || !f.exists() || !f.canRead()) return;
    	if(f.getName().startsWith(".") && f.getName().length()==1) return;
    	LOG.info(f);
    	if(f.isDirectory() && f.canRead())
    		{
    		if(f.getName().toLowerCase().equals("tmp")) return;
    		if(f.getName().toLowerCase().equals("jeter")) return;
    		
    		
    		if(f.getName().equals("Intensities") && f.getParentFile()!=null &&
    				f.getParentFile().getName().equals("Data"))
    			{
    			LOG.info("Skipping "+f);
    			return;
    			}
    		if(f.getName().startsWith("L") && f.getParentFile()!=null &&
    				(f.getParentFile().getName().equals("Thumbnail_Images") || 
    				 f.getParentFile().getName().equals("Processed")
    				))
    			{
    			LOG.info("Skipping "+f);
    			return;
    			}
    		
    		long now=System.currentTimeMillis();
    		if(now-lastDirTimeMillis > 30*1000)
	    		{
    			LOG.info("In "+f);
	    		this.lastDirTimeMillis=now;
	    		}
    		
    		
    		final File array[]=f.listFiles(this.fileFilter);
    		if(array!=null)
    			{
    			for(File f2:array) recursive(f2);    		
    			}
    		
    		
    		final Counter<String> fastqSamples=getFastqSampleInDirectory(f);
    		if(!fastqSamples.isEmpty())
    			{
    			fastqDir(f, fastqSamples);
    			}
    		}
    	else if(f.isFile() && this.fileFilter.accept(f))
    		{
    		analyze(f);
    		}
    	}
    
    private boolean isEntryShouldBeDeleted(File f)
    	{
    	if(!f.exists()) return true;
    	if(f.isDirectory() && getFastqSampleInDirectory(f).isEmpty()) return true;
    	return false;
    	}
    

    
    @Parameter(names={"-B","--bdb-home"},description="berkeleydb home directory",required=true)
    private File bdbHome=null;
    @Parameter(names={"-D","--dump"},description="dump as XML to stdout and exit")
    private boolean dump=false;
    
    @Override
    public int doWork(final List<String> args) {
		EnvironmentConfig envCfg=new EnvironmentConfig();
		Environment env=null;
		
		
		
		if(bdbHome==null)
			{
			LOG.error("BDB home undefined");
			return -1;
			}
		if(!bdbHome.exists() || !bdbHome.isDirectory()) {
			LOG.error("BDB doesn't exist or is not a directory");
			return -1;
			}
		
		File root=new File("/");
		Cursor cursor=null;
		OutputStream xmlout=null;
		try
			{
			if(dump)
				{
				if(!args.isEmpty()) {
					LOG.error("dump/illegal number of arguments");
					return -1;
					}
				}
			else {
			     if(args.size()==1)
				{
				root=new File(args.get(0));
				}
			      else
				{
				LOG.error("expected one file for root");
				return -1;
				}
			}
			
			envCfg.setAllowCreate(!dump);
			envCfg.setReadOnly(dump);
			envCfg.setTransactional(false);
			
			
			LOG.info("Opening env "+bdbHome);
			env=new Environment(bdbHome, envCfg);
			
			
			//TransactionConfig txnCfg=new TransactionConfig();
			//this.txn=env.beginTransaction(null, txnCfg);

			
			LOG.info("Opening database "+DATABASE_NAME);
			DatabaseConfig cfg=new DatabaseConfig();
			cfg.setAllowCreate(!dump);
			cfg.setReadOnly(dump);
			cfg.setTransactional(false);
			this.database=env.openDatabase(this.txn,DATABASE_NAME, cfg);
			
			DatabaseEntry key=new DatabaseEntry();
			DatabaseEntry data=new DatabaseEntry();

			if(dump)
				{
				XMLEventFactory xef=XMLEventFactory.newFactory();
	    		XMLOutputFactory xof=XMLOutputFactory.newFactory();
	    		XMLInputFactory xif=XMLInputFactory.newFactory();
	    		xmlout  = openFileOrStdoutAsStream(outputFile);
	    		XMLEventWriter out=xof.createXMLEventWriter(xmlout,"UTF-8");
	    		
	    		out.add(xef.createStartDocument("UTF-8", "1.0"));
	    		out.add(xef.createStartElement("", "","ngs-files"));
	    		cursor=this.database.openCursor(this.txn, null);
	    		while(cursor.getNext(key, data, LockMode.DEFAULT)==OperationStatus.SUCCESS)
	    			{
	    			final File file=new File(StringBinding.entryToString(key));
	    			if(isEntryShouldBeDeleted(file))
	    				{
	    				LOG.info("deleting entry for "+file);
	    				//cursor.delete();//no env is read only
	    				continue;
	    				}
	    			final StringReader sr=new StringReader(StringBinding.entryToString(data));
	    			XMLEventReader xr=xif.createXMLEventReader(sr);
	    			while(xr.hasNext())
	    				{
	    				XMLEvent evt=xr.nextEvent();
	    				if(evt.isStartDocument() ) continue;
	    				if(evt.isEndDocument()) break;
	    				out.add(evt);
	    				}
	    			xr.close();
	    			sr.close();
	    			}
	    		cursor.close();
	    		
	    		out.add(xef.createEndElement("","","ngs-files"));
	    		out.add(xef.createEndDocument());
	    		out.flush();
	    		out.close();
	    		xmlout.flush();
	    		xmlout.close();
	    		xmlout=null;
				}
			else
				{
				recursive(root);
				
				//final cleanup
				cursor=this.database.openCursor(this.txn, null);
				while(cursor.getNext(key, data, LockMode.DEFAULT)==OperationStatus.SUCCESS)
	    			{
	    			File file=new File(StringBinding.entryToString(key));
	    			
	    			if(isEntryShouldBeDeleted(file))
	    				{
	    				LOG.info("deleting entry for "+file);
	    				cursor.delete();
	    				}
	    			}
				cursor.close();
				}
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			
			if(cursor!=null)  try { cursor.close();} catch(Exception err){}
			if(this.txn!=null)try { this.txn.commit();} catch(Exception err){}
			if(this.database!=null) try { this.database.close();} catch(Exception err){}
			if(env!=null) try { env.close();} catch(Exception err){}
			CloserUtil.close(xmlout);
			}
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		/* https://bugs.openjdk.java.net/browse/JDK-8028111 */
		System.setProperty("jdk.xml.entityExpansionLimit","0");
		new NgsFilesScanner().instanceMainWithExit(args);

	}

}
