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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.PrintWriter;
import java.util.Collection;
import java.util.regex.Pattern;

import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.LineReader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

/**
 * Extends a BED by 'X' bases. Deprecated: use bedtools slop
 * @author lindenb
 *
 */
/** use bedtools slop */
@Deprecated //"use bedtools slop"
public class ExtendBed extends AbstractExtendBed
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(ExtendBed.class);


	 @Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  class MyCommand extends AbstractExtendBed.AbstractExtendBedCommand
	 	{		
		private SAMSequenceDictionary samSequenceDictionary = null;
	
    
		@Override
		protected Collection<Throwable> call(String inputName) throws Exception
				{
				if(REF==null)
					{
					return wrapException("undefined ref");
					}
				if(EXTEND<0)
					{
					return wrapException("negative value for X");
					}
				LineReader r=null;
				PrintWriter pw= null;
				try {
					
					this.samSequenceDictionary=new SAMSequenceDictionaryFactory().load(REF);
					if(inputName==null)
						{
						LOG.info("reading from stdin");
						r=new AsciiLineReader(stdin());
						}
					else
						{
						r=new AsciiLineReader(IOUtils.openURIForReading(inputName));
						}
					pw = super.openFileOrStdoutAsPrintWriter();
					
					final Pattern tab=Pattern.compile("[\t]");
					String line;
					while((line=r.readLine())!=null)
						{
						String tokens[]=tab.split(line);
						if(tokens.length<3)
							{
							return wrapException("Not enough cols in "+line);
							}
						SAMSequenceRecord rec;
						if((rec=this.samSequenceDictionary.getSequence(tokens[0]))==null)
							{
							return wrapException("Chromosome "+tokens[0]+" not in REF dictionary. ignoring");
							}
						int start=Integer.parseInt(tokens[1]);
						if(start<0)
							{
							return wrapException("start<0 in "+line );
							}
						int end=Integer.parseInt(tokens[2]);
						if(start>end)
							{
							return wrapException("  end<start in "+line );
							}
						start=Math.max(0, start-EXTEND);
						end=Math.min(end+EXTEND,rec.getSequenceLength());
						for(int i=0;i< tokens.length;++i)
							{
							if(i>0) pw.print('\t');
							switch(i)
									{
									case 1: pw.print(start); break;
									case 2: pw.print(end); break;
									default: pw.print(tokens[i]); break;
									}
							}
						pw.println();
						}	
					return RETURN_OK;
					} 
				catch (Exception e)
					{
					return wrapException(e);
					}
				finally
					{
					if(r!=null) r.close();
					pw.flush();
					pw.close();
					this.samSequenceDictionary=null;
					}
			}
	 	}
	
	public static void main(String[] args) {
		new ExtendBed().instanceMainWithExit(args);
		}
	}
