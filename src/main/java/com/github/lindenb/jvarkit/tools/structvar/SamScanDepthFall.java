/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.math.stats.Percentile;
import com.github.lindenb.jvarkit.tools.misc.ConcatSam;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.iterator.FilterIterator;
import com.github.lindenb.jvarkit.util.iterator.ForwardPeekIterator;
import com.github.lindenb.jvarkit.util.iterator.ForwardPeekIteratorImpl;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
/**

BEGIN_DOC

END_DOC
*/
@Program(name="samscandepthfall",
	description="Explore balanced translocations between two chromosomes using discordant paired-end reads.",
	keywords={"sam","bam","xslt","xml"}
	)
public class SamScanDepthFall extends Launcher {
	private static final Logger LOG = Logger.build(SamScanDepthFall.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Override
	public int doWork(final List<String> args) {
		final List<SamReader> samReaders = new ArrayList<>();
		try {
			final List<File> filenames = IOUtils.unrollFiles2018(args);
			SAMSequenceDictionary dict=null;
			for(int tid=0;tid< dict.size();++tid)
				{
				final SAMSequenceRecord ssr= dict.getSequence(tid);
				int depth[] = new int[ssr.getSequenceLength()];
				for(final SamReader reader: samReaders)
					{
					Arrays.fill(depth,0);
					final SAMRecordIterator iter = reader.query(ssr.getSequenceName(), 1,ssr.getSequenceLength(),false);
					while(iter.hasNext())
						{
						final SAMRecord rec= iter.next();
						if(rec.getReadUnmappedFlag()) continue;
						final Cigar cigar = rec.getCigar();
						if(cigar==null || cigar.isEmpty()) continue;
						int ref=rec.getStart();
						for(final CigarElement ce: cigar.getCigarElements())
							{
							final CigarOperator op=ce.getOperator();
							if(op.consumesReferenceBases())
								{
								if(op.consumesReadBases())
									{
									for(int i=0;i< ce.getLength();++i)
										{
										int pos=ref+i-1;
										if(pos >=0 && pos< depth.length)
											{
											depth[pos]++;
											}
										}
									}
								ref+=ce.getLength();
								}
							}
						}
					iter.close();
					}
				System.gc();
				}
			
			
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(samReaders);
			}
		}
	
	public static void main(String[] args) {
		new SamScanDepthFall().instanceMainWithExit(args);

	}

}
