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


*/
package com.github.lindenb.jvarkit.tools.biostar;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloserUtil;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamRecordFilterFactory;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;

/*
BEGIN_DOC
## Example

```bash
$ java -jar dist/sortsamrefname.jar --samoutputformat BAM input.bam  |\
  java -jar dist/biostar154220.jar  -n 20 --samoutputformat BAM |\
  samtools sort -T tmp -o output.bam -


$ samtools mpileup output.bam  | cut -f 4 | sort | uniq -c

  12692 0
 596893 1
  94956 10
  56715 11
  76947 12
  57912 13
  66585 14
  51961 15
  63184 16
  47360 17
  65189 18
  65014 19
 364524 2
 169064 20
  72078 3
 118288 4
  54802 5
  82555 6
  53175 7
  78474 8
  54052 9

```

END_DOC

*/
@Program(name="biostar154220",
	description="Cap BAM to a given coverage",
	biostars=154220,
	keywords={"bam","sam","coverage"},
	modificationDate="20190417"
	)
public class Biostar154220 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar154220.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-n","--depth"},description="number of reads")
	private int capDepth=20;
	
	@Parameter(names={"-filter","--filter"},description=SamRecordFilterFactory.FILTER_DESCRIPTION,converter=SamRecordFilterFactory.class,splitter=NoSplitter.class)
	private SamRecordFilter filter  = SamRecordFilterFactory.getDefault();

	@ParametersDelegate
	private WritingBamArgs writingBams=new WritingBamArgs();
	
	private int doWork(final SamReader in) throws IOException
		{
		final SAMFileHeader header= in.getFileHeader();
		if(header.getSortOrder()!=SAMFileHeader.SortOrder.unsorted)
			{
			LOG.error("input should be unsorted, reads sorted on REF/query-name e.g: see https://github.com/lindenb/jvarkit/wiki/SortSamRefName");
			return -1;
			}
		final SAMSequenceDictionary dict=SequenceDictionaryUtils.extractRequired(header);
		SAMFileWriter out=null;
		SAMRecordIterator iter=null;
		int prev_tid=-1;
		int depth_array[]=null;
		try
			{
			SAMFileHeader header2=header.clone();
			header2.addComment("Biostar154220"+" "+getVersion()+" "+getProgramCommandLine());
			out = this.writingBams.openSAMFileWriter(outputFile,header2, true);
			final ProgressFactory.Watcher<SAMRecord> progress=ProgressFactory.newInstance().dictionary(dict).logger(LOG).build();
			iter = in.iterator();
			List<SAMRecord> buffer=new ArrayList<>();
			for(;;)
				{
				SAMRecord rec =null;
				
				if(iter.hasNext())
					{
					rec = progress.apply(iter.next());
					}
				
				if(rec!=null && rec.getReadUnmappedFlag())
					{
					out.addAlignment(rec);
					continue;
					}
				//no more record or 
				if(!buffer.isEmpty() &&
					rec!=null &&
					buffer.get(0).getReadName().equals(rec.getReadName()) &&
					buffer.get(0).getReferenceIndex().equals(rec.getReferenceIndex())
					)
					{
					buffer.add(progress.apply(rec));
					}
				else if(buffer.isEmpty() && rec!=null)
					{
					buffer.add(progress.apply(rec));
					}
				else //dump buffer
					{
					if(!buffer.isEmpty())
						{
						final int tid = buffer.get(0).getReferenceIndex();
						if(prev_tid==-1 || prev_tid!=tid)
							{
							SAMSequenceRecord ssr=dict.getSequence(tid);
							prev_tid=tid;
							depth_array=null;
							System.gc();
							LOG.info("Alloc memory for contig "+ssr.getSequenceName()+" N="+ssr.getSequenceLength()+"*sizeof(int)");
							depth_array=new int[ssr.getSequenceLength()+1];//use a +1 pos
							Arrays.fill(depth_array, 0);
							}
						//position->coverage for this set of reads
						Counter<Integer> readposition2coverage=new Counter<Integer>();
						
						boolean dump_this_buffer=true;
						for(final SAMRecord sr:buffer)
							{
							if(!dump_this_buffer) break;
							if(this.filter.filterOut(sr)) continue;
							
							
							final Cigar cigar=sr.getCigar();
							if(cigar==null)
								{
								throw new IOException("Cigar missing in "+rec.getSAMString());
								}
							int refPos1=sr.getAlignmentStart();
							for(final CigarElement ce:cigar.getCigarElements())
								{
								final CigarOperator op =ce.getOperator();
								if(!op.consumesReferenceBases()) continue;
								if(op.consumesReadBases())
									{
									for(int x=0;x<ce.getLength() && refPos1+x< depth_array.length;++x)
										{
										int cov = (int)readposition2coverage.incr(refPos1+x);
										if( depth_array[refPos1+x]+cov > this.capDepth)
											{
											dump_this_buffer=false;
											break;
											}
										}
									}
								if(!dump_this_buffer) break;
								refPos1+=ce.getLength();
								}
							}
						if(dump_this_buffer)
							{
							//consumme this coverage
							for(Integer pos:readposition2coverage.keySet())
								{
								depth_array[pos]+= (int)readposition2coverage.count(pos);
								}
							for(SAMRecord sr:buffer)
								{
								out.addAlignment(sr);
								}
							}
						
						buffer.clear();
						}
					if(rec==null) break;
					buffer.add(rec);
					}
				}
			depth_array=null;
			progress.close();
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(out);
			}
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		if(this.capDepth<0) // -1 == infinite
			{
			LOG.error("Bad depth:"+this.capDepth);
			return -1;
			}
		SamReader in=null;
		try
			{
			in=openSamReader(oneFileOrNull(args));
			return doWork(in); 
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			}
		}

	
	public static void main(final String[] args) throws IOException
		{
		new Biostar154220().instanceMainWithExit(args);
		}
		

	}
