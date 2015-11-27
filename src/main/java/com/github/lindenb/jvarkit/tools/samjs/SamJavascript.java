/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.samjs;



/**
 * Author: Pierre Lindenbaum PhD. @yokofakun
 * Motivation http://www.biostars.org/p/66319/ 
 */

import java.util.Collection;

import javax.script.Bindings;
import javax.script.CompiledScript;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;



public class SamJavascript
	extends AbstractSamJavascript
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(SamJavascript.class);

	private CompiledScript  script=null;
	private SAMFileWriter failingReadsWriter=null;
	

	public SamJavascript()
		{
		
		}
	
	/* open failing bam if it was not already open */
	private void openFailing(SAMFileHeader h)
		{
		if(this.failingReadsFile==null) return;
		if(this.failingReadsWriter==null)
			{
			LOG.info("Writing failings to "+ this.failingReadsFile);
			SAMFileHeader h2= h.clone();
			SAMFileWriterFactory sfwf = new SAMFileWriterFactory();
			sfwf.setCreateMd5File(false);
			sfwf.setCreateIndex(false);
		
			this.failingReadsWriter=sfwf.makeSAMOrBAMWriter(
					h2, 
					true
					,this.failingReadsFile);
			}
		}
	
	private void failing(SAMRecord rec,SAMFileHeader h)
		{
		openFailing(h);
		if(failingReadsWriter!=null) failingReadsWriter.addAlignment(rec);
		}
	
	

	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		
		 SAMRecordIterator iter=null;
		SamReader samFileReader=null;
		SAMFileWriter sw=null;
		try
			{
			this.script  = super.compileJavascript();
			samFileReader= openSamReader(inputName);
			final SAMFileHeader header=samFileReader.getFileHeader();
			sw = openSAMFileWriter(header, true);
			
			
			long count=0L;
	        Bindings bindings = this.script.getEngine().createBindings();
	        bindings.put("header", samFileReader.getFileHeader());
	        SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
	        iter = samFileReader.iterator();
			while(iter.hasNext())
				{
				SAMRecord record=iter.next();
				progress.watch(record);
				bindings.put("record", record);
				
				if(super.evalJavaScriptBoolean(this.script, bindings))
					{
					++count;
					sw.addAlignment(record);
					if(this.LIMIT>0L && count>=this.LIMIT) break;
					
					}
				else
					{
					failing(record,header);
					}
				}
			sw.close();
			/* create empty if never called */
			openFailing(header);
			
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(samFileReader);
			CloserUtil.close(sw);
			CloserUtil.close(failingReadsWriter);
			}
		}
	
		
	public static void main(String[] args) throws Exception
		{
		new SamJavascript().instanceMainWithExit(args);
		}
	
	}
