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
package com.github.lindenb.jvarkit.tools.tview;

import java.io.File;
import java.io.PrintWriter;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.IntervalUtils;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;


public class TViewCmd extends AbstractCommandLineProgram
	{
	private File referenceFile=null;
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	
	private abstract class Node
		{
		Node parent=null;
		Node next=null;
		Node child=null;
		public void appendChild(Node n)
			{
			if(child==null)
				{
				this.child=n;
				}
			else
				{
				Node curr=child;
				while(curr.next!=null) curr=curr.next;
				curr.next=n;
				}
			n.parent=this;
			}

		}
	private class Element
		extends Node
		{
		}
	
	private class Root
		extends Node
		{
		Root()
			{
			}
		}
	
	private void scan()
		{
		String chrom="";
		int chromStart=0;
		int width=80;
		int column2genomic[]=new int[width];
		for(int i=0;i< column2genomic.length;++i)
			{
			column2genomic[i] = chromStart+i;
			}
		GenomicSequence genomicSequence=new GenomicSequence(indexedFastaSequenceFile, chrom);
		Element contig=new Element();
		for(int i=0;i< width;++i)
			{
			
			}
		
		Root root=new Root();
		SAMRecordIterator iter=null;
		while(iter.hasNext())
			{
			SAMRecord rec = iter.next();
			}
		}
	
	
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"Tview";
		}
    
    @Override
	public String getProgramDescription() {
		return "java equivalent of samtools tview";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println("-p (chrom:pos) region");
		super.printOptions(out);
		}
	
	public void setRegion(String chrom,int start)
		{
		
		}
	
	@Override
	public int doWork(String[] args)
		{
		String region=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args, getGetOptDefault()+"r:R:"))!=-1)
			{
			switch(c)
				{
				case 'r': region=opt.getOptArg();break;
				case 'R': this.referenceFile = new File(opt.getOptArg());break;
				default: 
					{
					switch(handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default: break;
						}
					}
				}
			}
		
		File faidx=null;
		File bamFile=null;
		
		if(opt.getOptInd()+1==args.length )
			{
			bamFile=new File(args[opt.getOptInd()]);
			}
		else if(opt.getOptInd()+2==args.length )
			{
			bamFile=new File(args[opt.getOptInd()]);
			faidx=new File(args[opt.getOptInd()+1]);
			}
		else
			{
			System.err.println("Illegal Number of arguments.");
			return -1;
			}
		
		
	
		
    	IndexedFastaSequenceFile ref=null;
		PrintWriter out=new PrintWriter(System.out);
		SamReader samReader=null;
		try {
			info("opening "+bamFile);
			SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
	        samReader=srf.open(bamFile);
	        
			Interval interval=IntervalUtils.parseOne(
					samReader.getFileHeader().getSequenceDictionary(),
					region
					);
			if(interval==null)
				{
				error("Bad interval "+interval);
				return -1;
				}
			if(faidx!=null)
				{
				ref=new IndexedFastaSequenceFile(faidx);
				}
	  
	        TViewHandler handler=new AsciiHandler();
	        new TView().execute(samReader, ref, interval, handler);
			} 
		catch (Exception e) {
			e.printStackTrace();
			return -1;
			}
		finally
			{
			out.flush();
			CloserUtil.close(samReader);
			CloserUtil.close(ref);
			}
		return 0;
		}

public static void main(String[] args)
	{
	new TViewCmd().instanceMainWithExit(args);
	}
}
