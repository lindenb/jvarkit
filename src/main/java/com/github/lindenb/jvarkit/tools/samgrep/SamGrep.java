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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.samgrep;


import java.io.BufferedReader;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;

public class SamGrep extends AbstractSamGrep
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(SamGrep.class);

	private final Map<String,Integer> readNames=new HashMap<String,Integer>(); 
    private SamGrep()
    	{
    	}
    
    @Override
    public Collection<Throwable> initializeKnime() {
    	readNames.clear();
    	
    	if(namefile!=null) {
	    	BufferedReader in=null;
			try
				{
				in=IOUtils.openFileForBufferedReading(super.namefile);
		    	String line;
		    	while((line=in.readLine())!=null)
		    		{
		    		line=line.trim();
		    		if(line.isEmpty()) continue;
		    		readNames.put(line,0);
		    		}
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(in);
				}
	    	}
    	for(final String line: super.nameStrings) {
    		readNames.put(line,0);
    		}
    	if(readNames.isEmpty())
			{
			LOG.warn("no read found.");
			}
    	return super.initializeKnime();
    }
    
    @Override
    public void disposeKnime() {
    	readNames.clear();
    	super.disposeKnime();
    }
    
	@Override
	protected Collection<Throwable> call(final String inputName) throws Exception {
		SAMFileWriter sfw=null;
		SAMFileWriter samStdout=null;
		SamReader sfr=null;
		try {
			sfr = super.openSamReader(inputName);
			final SAMFileHeader header=sfr.getFileHeader().clone();
			final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header);

			if(super.outputFile==null)
				{
				sfw= super.openSAMFileWriter(header, true);
				}
			else
				{
				if(super.divertToStdout) {
					if(SamReader.Type.SAM_TYPE.equals(super.samoutputformat)) {
					samStdout= super.createSAMFileWriterFactory().makeSAMWriter(header, true,stdout());
					} else
					{
					samStdout= super.createSAMFileWriterFactory().makeBAMWriter(header, true,stdout());
					}
				}
				sfw= super.openSAMFileWriter(header, true);
				}
		
			SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				boolean keep=false;
				final SAMRecord rec=progress.watch(iter.next());
				if(samStdout!=null) samStdout.addAlignment(rec);
				Integer count = readNames.get(rec.getReadName());
				if(count!=null)
					{
					keep=true;
					}
				if(super.inverse) keep=!keep;
				if(keep)
					{
					sfw.addAlignment(rec);
					}
				
				if(n_before_remove!=-1 && !inverse && keep)
					{
					count++;
					if(count>=n_before_remove)
						{
						readNames.remove(rec.getReadName());
						if(samStdout==null && readNames.isEmpty()) break;
						}
					else
						{
						readNames.put(rec.getReadName(),count);
						}
					}
				}
			progress.finish();
			return RETURN_OK;
		} catch (final Exception e) {
			return wrapException(e);
		} finally {
			CloserUtil.close(samStdout);
			CloserUtil.close(sfw);
			CloserUtil.close(sfr);
			}
		}
    
    public static void main(final String[] argv)
		{
	    new SamGrep().instanceMainWithExit(argv);
		}	

	}
