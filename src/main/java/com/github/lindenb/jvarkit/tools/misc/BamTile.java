
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
* 2015 : creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;

public class BamTile
	extends AbstractBamTile
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(BamTile.class);
	public BamTile()
			{
			}
	

	 @Override
		public  Command createCommand() {
				return new MyCommand();
			}
			 
		private static class MyCommand extends AbstractBamTile.AbstractBamTileCommand
			 	{
				@Override
				public Collection<Throwable> call() throws Exception 
					{
					List<String> args= this.getInputFiles();	
							
					SAMRecordIterator iter=null;
					SamReader sfr=null;
					SAMFileWriter sfw=null;
					try
						{
						
						SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
						
						if(args.isEmpty())
							{
							LOG.info("Reading from stdin");
							sfr=srf.open(SamInputResource.of(stdin()));
							}
						else if(args.size()==1)
							{
							File fin=new File(args.get(0));
							LOG.info("Reading from "+fin);
							sfr=SamFileReaderFactory.mewInstance().open(fin);
							}
						else
							{
							return wrapException("Illegal number of arguments.");
							}
						
						SAMFileHeader header1=sfr.getFileHeader();
						if(header1==null)
							{
							return wrapException("File header missing");
							}
						
						if(header1.getSortOrder()!=SAMFileHeader.SortOrder.coordinate)
							{
							return wrapException("File header not sorted on coordinate");
							}
						
						SAMFileHeader header2=header1.clone();
						header2.addComment(getName()+":"+getVersion()+":"+getProgramCommandLine());
							
						SAMFileWriterFactory sfwf=new SAMFileWriterFactory();
						sfwf.setCreateIndex(false);
						sfwf.setCreateMd5File(false);
						if(getOutputFile()==null)
							{
							sfw = sfwf.makeSAMWriter(header2, true, stdout());
							}
						else
							{
							sfw = sfwf.makeSAMOrBAMWriter(header2, true, getOutputFile());
							}
						
						
						
						SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header1);
						iter=sfr.iterator();
						LinkedList<SAMRecord> buffer=new LinkedList<>();
						for(;;)
							{
							SAMRecord rec=null;
							if( iter.hasNext())
								{
								rec= progress.watch(iter.next());
								if(rec.getReadUnmappedFlag()) continue;
								if(!buffer.isEmpty())
									{
									SAMRecord last= buffer.getLast();
									if( last.getReferenceIndex()==rec.getReferenceIndex() &&
										last.getAlignmentStart()<=rec.getAlignmentStart() &&
										last.getAlignmentEnd()>=rec.getAlignmentEnd())
										{
										continue;
										}
									}
								}
							if(rec==null || (!buffer.isEmpty() && buffer.getLast().getReferenceIndex()!=rec.getReferenceIndex()))
								{
								while(!buffer.isEmpty())
									{
									sfw.addAlignment(buffer.removeFirst());
									}
								if(rec==null) break;
								}
							buffer.add(rec);
							
							
							if(buffer.size()>2)
								{
								int index = buffer.size();
								SAMRecord prev =  buffer.get(index-3);
								SAMRecord curr =  buffer.get(index-2);
								SAMRecord next =  buffer.get(index-1);
								
								if( prev.getAlignmentEnd() >= next.getAlignmentStart() ||
									curr.getAlignmentEnd() <= prev.getAlignmentEnd())
									{
									buffer.remove(index-2);
									}
								else if(curr.getAlignmentStart() == prev.getAlignmentStart() &&
										curr.getAlignmentEnd() > prev.getAlignmentEnd()
										)
									{
									buffer.remove(index-3);
									}
								
								}
							while(buffer.size()>3)
								{
								sfw.addAlignment(buffer.removeFirst());
								}
							
							}
						progress.finish();
						return Collections.emptyList();
						}
					catch(Exception err)
						{
						return wrapException(err);
						}
					finally
						{
						CloserUtil.close(iter);
						CloserUtil.close(sfr);
						CloserUtil.close(sfw);
						}
					}
				}
	

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new BamTile().instanceMainWithExit(args);
		}
	}
