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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.OtherCanonicalAlign;
import com.github.lindenb.jvarkit.util.picard.OtherCanonicalAlignFactory;

public class SamShortInvertion extends AbstractSamShortInvertion
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(AbstractSamShortInvertion.class);

	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		SamReader r=null;
		PrintStream out = null;
		//SAMFileWriter w=null;
		try
			{
			r = openSamReader(inputName);
			out =  openFileOrStdoutAsPrintStream();
			final SAMFileHeader header=r.getFileHeader();
			OtherCanonicalAlignFactory xpalignFactory=new OtherCanonicalAlignFactory(header);
			int prev_tid=-1;
			short coverage[]=null;
			short max_coverage=0;
			//w=swf.make(header, System.out);
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			final SAMRecordIterator it= r.iterator();
			for(;;)
				{
				SAMRecord rec=null;
				if(it.hasNext())
					{
					rec=progress.watch(it.next());
					if(rec.getReadUnmappedFlag()) continue;
					if(rec.isSecondaryOrSupplementary()) continue;
					if(rec.getDuplicateReadFlag()) continue;
					if(rec.getReadFailsVendorQualityCheckFlag()) continue;
					}
				
				if(rec==null || prev_tid==-1 || (prev_tid!=-1 && prev_tid!=rec.getReferenceIndex()))
					{
					if(coverage!=null)
						{

						while(max_coverage>=Math.max(1,this.min_coverage))
							{
							LOG.info("Scanning "+header.getSequence(prev_tid).getSequenceName()+" for cov:"+max_coverage);

							int chromStart0=0;
							while(chromStart0 < coverage.length)
								{
								if(coverage[chromStart0]==max_coverage)
									{
									coverage[chromStart0]=0;//reset that pos
									int chromEnd0=chromStart0+1;
									while(chromEnd0 < coverage.length && coverage[chromEnd0]==max_coverage)
										{
										coverage[chromEnd0]=0;//reset that pos
										++chromEnd0;
										}
									out.print(header.getSequence(prev_tid).getSequenceName());
									out.print('\t');
									out.print(chromStart0);
									out.print('\t');
									out.print(chromEnd0);
									out.print('\t');
									out.print(max_coverage);
									out.println();
									
									//reset 3'
									for(int x=chromEnd0;x<coverage.length && coverage[x]>0;++x)
										{
										coverage[x]=0;
										}
									//reset 5'
									for(int x=chromStart0-1;x>=0 && coverage[x]>0;--x)
										{
										coverage[x]=0;
										}
									chromStart0=chromEnd0;
									}
								else
									{
									++chromStart0;
									}
								}
							max_coverage--;
							}
						coverage=null;
						}
					
					if(rec==null) break;
					prev_tid=rec.getReferenceIndex();
					LOG.info("Alloc sizeof(short)*"+header.getSequence(prev_tid).getSequenceLength());
					coverage=new short[header.getSequence(prev_tid).getSequenceLength()];
					Arrays.fill(coverage,(short)0);
					max_coverage=0;
					}
				
				List<OtherCanonicalAlign> saList=xpalignFactory.getXPAligns(rec);
				if(saList.isEmpty()) continue;
				for(OtherCanonicalAlign xp:saList)
					{
					if(!xp.getReferenceName().equals(rec.getReferenceName())) continue;
					
					if(!rec.getReadNegativeStrandFlag()) //read is plus
						{
						if(!xp.getReadNegativeStrandFlag())
							{
							//info(xp+" vs strand "+rec.getReadNegativeStrandFlag());
							continue;//ignore both same strand
							}
						}
					else //read.strand=='-'
						{
						if(xp.getReadNegativeStrandFlag())
							{
							//info(xp+" vs strand "+rec.getReadNegativeStrandFlag());
							continue;//ignore both same strand
							}
						}
					if(Math.abs(rec.getUnclippedStart() - xp.getAlignmentStart()) > max_size_inversion) 
						{
						//info(xp+" vs pos "+rec.getUnclippedStart());
						continue;
						}
					int chromStart=Math.min(rec.getUnclippedStart(), xp.getAlignmentStart())-1;
					int chromEnd=Math.max(rec.getUnclippedEnd(), xp.getAlignmentStart())-1;
					for(int x=chromStart;x<=chromEnd && x< coverage.length;++x )
						{
						if(coverage[x]<Short.MAX_VALUE) coverage[x]++;
						if(max_coverage< coverage[x])
							{
							LOG.info("Max coverage "+max_coverage);
							max_coverage=coverage[x];
							}
						}
					}
				}
			
			it.close();
			progress.finish();
			return RETURN_OK;
			}

		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(out);
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new SamShortInvertion().instanceMainWithExit(args);
		}
	}
