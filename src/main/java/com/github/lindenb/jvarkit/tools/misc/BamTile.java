
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
import com.github.lindenb.jvarkit.knime.AbstractKnimeApplication;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;

public class BamTile
	extends AbstractKnimeApplication
	{
	public BamTile()
			{
			}
	
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"BamTile";
		}
	
	@Override
	public String getProgramDescription() {
		return "Answer to @sjackman : Is there a bedtools command to determine a minimal tiling path? A minimal set of features that cover a maximum of the target. ";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -o (output file) default: SAM as stdout");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:"))!=-1)
			{
			switch(c)
				{
				case 'o': this.setOutputFile(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		return mainWork(opt.getOptInd(), args);
		}
	
	@Override
	public int executeKnime(List<String> args)
		{
		SAMRecordIterator iter=null;
		SamReader sfr=null;
		SAMFileWriter sfw=null;
		try
			{
			SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
			
			if(args.isEmpty())
				{
				info("Reading from stdin");
				sfr=srf.open(SamInputResource.of(System.in));
				}
			else if(args.size()==1)
				{
				File fin=new File(args.get(0));
				info("Reading from "+fin);
				sfr=SamFileReaderFactory.mewInstance().open(fin);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			
			SAMFileHeader header1=sfr.getFileHeader();
			if(header1==null)
				{
				error("File header missing");
				return -1;
				}
			
			if(header1.getSortOrder()!=SAMFileHeader.SortOrder.coordinate)
				{
				error("File header not sorted on coordinate");
				return -1;
				}
			
			SAMFileHeader header2=header1.clone();
			header2.addComment(getProgramName()+":"+getVersion()+":"+getProgramCommandLine());
				
			SAMFileWriterFactory sfwf=new SAMFileWriterFactory();
			sfwf.setCreateIndex(false);
			sfwf.setCreateMd5File(false);
			if(getOutputFile()==null)
				{
				sfw = sfwf.makeSAMWriter(header2, true, System.out);
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
					
					}
				while(buffer.size()>3)
					{
					sfw.addAlignment(buffer.removeFirst());
					}
				
				}
			progress.finish();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(sfr);
			CloserUtil.close(sfw);
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
