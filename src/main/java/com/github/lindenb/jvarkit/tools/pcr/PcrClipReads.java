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


History:
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.pcr;

import java.io.BufferedReader;
import java.util.Collection;
import java.util.Iterator;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class PcrClipReads extends AbstractPcrClipReads
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(PcrClipReads.class);

	private IntervalTreeMap<Interval> bedIntervals=new IntervalTreeMap<Interval>();
	
	private Interval findInterval(final SAMRecord rec)
		{
		if(rec.getReadUnmappedFlag()) return null;
		return findInterval(rec.getContig(), rec.getAlignmentStart(), rec.getAlignmentEnd());
		}
	private Interval findInterval(String chrom,int start,int end)
		{
		Interval i= new Interval(chrom,start,end);
		Collection<Interval> L=this.bedIntervals.getOverlapping(i);
		Iterator<Interval> iter = L.iterator();
		if(iter.hasNext())
			{
			Interval j = iter.next();
			if(iter.hasNext()  ) throw new IllegalStateException("Overlapping PCR intervals : "+j+" "+iter.next());
			return j;
			}
		return null;
		}
	
	
	private Collection<Throwable> run(SamReader reader)
		{
		SAMFileHeader header1= reader.getFileHeader();
		SAMFileHeader header2 = header1.clone();
		header2.addComment(getName()+" "+getVersion()+": Processed with "+getProgramCommandLine());
		header2.setSortOrder(SortOrder.unsorted);
		SAMFileWriter sw=null;
		SAMRecordIterator iter = null;
		try
			{
			sw = super.openSAMFileWriter(header2, false);
			
			final SAMSequenceDictionaryProgress progress =new SAMSequenceDictionaryProgress(header1);
			iter =  reader.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec= progress.watch(iter.next());
				
				if(super.onlyFlag!=-1 &&  (rec.getFlags() & super.onlyFlag) != 0) {
					sw.addAlignment(rec);
					continue;
				}
				
				if(rec.getReadUnmappedFlag())
					{
					sw.addAlignment(rec);
					continue;
					}
				Interval fragment = findInterval(rec);
				if(fragment==null)
					{
					rec.setMappingQuality(0);
					sw.addAlignment(rec);
					continue;
					}
				// strand is '-' and overap in 5' of PCR fragment
				if( rec.getReadNegativeStrandFlag() &&
					fragment.getStart()< rec.getAlignmentStart() &&
					rec.getAlignmentStart()< fragment.getEnd())
					{
					rec.setMappingQuality(0);
					sw.addAlignment(rec);
					continue;
					}
				// strand is '+' and overap in 3' of PCR fragment
				if( !rec.getReadNegativeStrandFlag() &&
					fragment.getStart()< rec.getAlignmentEnd() &&
					rec.getAlignmentEnd()< fragment.getEnd())
					{
					rec.setMappingQuality(0);
					sw.addAlignment(rec);
					
					continue;
					}
				
				// contained int PCR fragment
				if(rec.getAlignmentStart()>= fragment.getStart() && rec.getAlignmentEnd()<=fragment.getEnd())
					{
					sw.addAlignment(rec);
					
					continue;
					}
				ReadClipper readClipper = new ReadClipper();
				rec = readClipper.clip(rec, fragment);		
				sw.addAlignment(rec);
				}
			progress.finish();
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(sw);
			}
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {

		if(super.bedFile==null)
			{
			return wrapException("undefined bed file option -"+OPTION_BEDFILE);
			}
		BufferedReader r=null;
		SamReader samReader=null;
		try {
			samReader = openSamReader(inputName);
			LOG.info("reading bed File "+super.bedFile);
			final Pattern tab= Pattern.compile("[\t]");
			r= IOUtils.openFileForBufferedReading(super.bedFile);
			String line;
			while((line=r.readLine())!=null)
				{
				String tokens[]=tab.split(line);
				if(tokens.length<3)
					{
					return wrapException("Bad bed line "+line);
					}
				String chrom = tokens[0];
				int chromStart1 = Integer.parseInt(tokens[1])+1;
				int chromEnd1 = Integer.parseInt(tokens[2])+0;
				if(chromStart1<1 || chromStart1>chromEnd1)
					{
					wrapException("Bad bed line "+line);
					}
				Interval i =new Interval(chrom, chromStart1, chromEnd1);
				this.bedIntervals.put(i, i);
				}
			
			CloserUtil.close(r);r=null;
			
			return run(samReader);
			}
		catch (Exception err) {
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(samReader);
			this.bedIntervals.clear();;
			}
		}

	
	public static void main(String[] args) {
		new PcrClipReads().instanceMain(args);
		}

}
