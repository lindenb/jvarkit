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
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamRecordFilterFactory;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloserUtil;


/**

BEGIN_DOC


Example:


```
$ curl -s "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20110915_CEUtrio_b37_decoy_alignment/CEUTrio.HiSeq.WGS.b37_decoy.NA12892.clean.dedup.recal.bam" |  java -jar dist/biostar234230.jar 
#contig	start	end	pairs_in_window	pairs_over_window	pairs_partial_overlap
1	10000	10100	0	2	240
1	10050	10150	4	615	274
1	10100	10200	0	800	276
1	10150	10250	0	216	649
1	10200	10300	0	2982	809
1	10250	10350	0	2918	207
1	10300	10400	0	1923	2851
1	10350	10450	0	227	4498
1	10400	10500	0	31	1971
(...)

```





END_DOC
*/


@Program(name="biostar234230",
	description="Sliding Window : discriminate partial and fully contained fragments (from a bam file)",
	biostars=234230,
	keywords= {"sam","bam"},
	modificationDate="20190417"
	)
public class Biostar234230 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar234230.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-w","--winsize"},description="Window size. "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int windowSize = 100 ;

	@Parameter(names={"-s","--winshift"},description="Shift each window by 's' bases." + DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int windowShift = 50 ;

	@Parameter(names={"-filter","--filter"},description=SamRecordFilterFactory.FILTER_DESCRIPTION,converter=SamRecordFilterFactory.class,splitter=NoSplitter.class)
	private SamRecordFilter filter  = SamRecordFilterFactory.getDefault();

	
	private static class SlidingWindow
		{
		final int start;
		final int end;
		int pairs_in_window=0;
		int pairs_over_window=0;
		int pairs_partial_overlap=0;
		
		SlidingWindow(int start,int end) {
		this.start = start;
		this.end=end;
		}
		
		void print(PrintWriter out ,final String contig) {
			out.print(contig);
			out.print("\t");
			out.print( start);
			out.print("\t");
			out.print( end);
			out.print("\t");
			out.print( pairs_in_window);
			out.print("\t");
			out.print( pairs_over_window);
			out.print("\t");
			out.print( pairs_partial_overlap);
			out.println();
			}
		
		void watch(int chromStart,int chromEnd) {
			if(chromStart>=start && chromEnd<=end) {
				++pairs_in_window;
				}
			else if(chromEnd<this.start || chromStart>this.end)
				{
				//nothing
				}
			else if(chromStart<this.start && chromEnd > this.end)
				{
				++pairs_over_window;
				}
			else
				{
				++pairs_partial_overlap;
				}
			}
		}
	
	@Override
	public int doWork(final List<String> args) {		
		
		if(this.windowSize<=0)
			{
			LOG.error("Bad window size.");
			return -1;
			}
		if(this.windowShift<=0)
			{
			LOG.error("Bad window shift.");
			return -1;
			}
		
		
		SamReader in = null;
		PrintWriter out=null;
		SAMRecordIterator iter=null;
		final List<SlidingWindow> buffer= new ArrayList<>();
		try
			{
			
			in = super.openSamReader(oneFileOrNull(args));
			final SAMFileHeader header = in.getFileHeader();
			if(header.getSortOrder()!=SAMFileHeader.SortOrder.coordinate) {
			LOG.error("Bam is not sorted on coordinate");
			return -1;
			}
			out = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			
			out.print("#contig");
			out.print("\t");
			out.print("start");
			out.print("\t");
			out.print("end");
			out.print("\t");
			out.print("pairs_in_window");
			out.print("\t");
			out.print("pairs_over_window");
			out.print("\t");
			out.print("pairs_partial_overlap");
			out.println();
			
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
			iter = in.iterator();
			String prev_contig=null;
			for(;;)
				{
				final SAMRecord rec = (iter.hasNext()?progress.watch(iter.next()):null);
				if(rec!=null )
					{
					if(rec.getReadUnmappedFlag()) continue;
					if(this.filter.filterOut(rec)) continue;
					if(rec.getReadUnmappedFlag()) continue;
					if(!rec.getReadPairedFlag()) continue;
					if(rec.getMateUnmappedFlag()) continue;
					if(!rec.getProperPairFlag()) continue;
					if(!rec.getReferenceName().equals(rec.getMateReferenceName())) continue;
					if(rec.getMateAlignmentStart()> rec.getAlignmentStart()) continue;
					
					}
				if(rec==null || !rec.getContig().equals(prev_contig))
					{
					for(final SlidingWindow w:buffer)
						{
						w.print(out,prev_contig);
						}
					buffer.clear();
					if(rec==null) break;
					prev_contig = rec.getContig();
					}
				final int fragStart = rec.getMateAlignmentStart();
				final int fragEnd = rec.getAlignmentEnd();
				while(!buffer.isEmpty() && buffer.get(0).end < fragStart)
					{
					buffer.remove(0).print(out,rec.getContig());;
					}
				final int winStart = fragStart - fragStart%this.windowSize;
				final int winEnd = fragEnd - fragEnd%this.windowSize;
				
				for(int winPos = winStart;winPos<=winEnd;winPos+=this.windowShift) {
					if(buffer.isEmpty() || buffer.get(buffer.size()-1).start < winPos)
						{
						buffer.add(new SlidingWindow(winPos,winPos+this.windowSize));
						}
					}
				for(SlidingWindow w:buffer)
					{
					w.watch(fragStart,fragEnd);
					}				
				}
			
				
			progress.finish();
			out.flush();
			return RETURN_OK;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(in);
			CloserUtil.close(out);
			}	
		}
	
	public static void main(final String[] args) {
		new Biostar234230().instanceMainWithExit(args);
		}
	}
