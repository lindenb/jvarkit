/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloseableIterator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassBamLauncher;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;

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

## Cited in

  * "Burden of Cardiomyopathic Genetic Variation in Lethal Pediatric Myocarditis". Amy R. Kontorovich & al. 6 Jul 2021 https://doi.org/10.1161/CIRCGEN.121.003426Circulation: Genomic and Precision Medicine. 2021;0


END_DOC

*/
@Program(name="biostar154220",
	description="Cap BAM to a given coverage",
	biostars={154220,9471266},
	keywords={"bam","sam","coverage","depth"},
	modificationDate="20210312",
	creationDate="20150812"
	)
public class Biostar154220 extends OnePassBamLauncher {
	private static final Logger LOG = Logger.build(Biostar154220.class).make();
	
	@Parameter(names={"-d","-n","--depth"},description="expected coverage.")
	private int capDepth=20;
	
	@Parameter(names={"-filter","--filter"},description=SamRecordJEXLFilter.FILTER_DESCRIPTION,converter=SamRecordJEXLFilter.StringConverter.class,splitter=NoSplitter.class)
	private SamRecordFilter filter  = SamRecordJEXLFilter.buildDefault();

	@Parameter(names={"--keep-unmapped"},description="write unmapped reads")
	private boolean keep_unmapped=false;
	@Parameter(names={"--query-sorted"},description="Input was sorted on query name but I promess there is one and only one chromosome: e.g: samtools view -h in.bam 'chr1:234-567' | samtools sort -n -) .")
	private boolean query_sorted=false;


	
	private SAMSequenceDictionary dict = null;
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int beforeSam()
		{
		if(this.capDepth<0) // -1 == infinite
			{
			LOG.error("Bad depth:"+this.capDepth);
			return -1;
			}
		return super.beforeSam();
		}
	
	@Override
	protected SAMFileHeader createOutputHeader(final SAMFileHeader headerIn) {
		this.dict=SequenceDictionaryUtils.extractRequired(headerIn);
		if(this.query_sorted) {
			if(headerIn.getSortOrder()!=SAMFileHeader.SortOrder.queryname) {
				throw new SAMException("input should be sorted on queryname, but got " + headerIn.getSortOrder());
				}
			}
		else if(headerIn.getSortOrder()!=SAMFileHeader.SortOrder.unsorted)
			{
			throw new SAMException("input should be unsorted, reads sorted on REF/query-name e.g: see https://github.com/lindenb/jvarkit/wiki/SortSamRefName");
			}
		final SAMFileHeader headerOut= super.createOutputHeader(headerIn);
		headerOut.setSortOrder(SAMFileHeader.SortOrder.unsorted);
		return headerOut;
		}
	
	@Override
	protected void scanIterator(final SAMFileHeader headerIn,
			CloseableIterator<SAMRecord> iter,
			final SAMFileWriter out)
		{
		int prev_tid=-1;
		int[] depth_array = null;
		final List<SAMRecord> buffer=new ArrayList<>();
		for(;;)
			{
			SAMRecord rec =null;
			
			if(iter.hasNext())
				{
				rec = iter.next();
				}
			
			if(rec!=null && rec.getReadUnmappedFlag())
				{
				if(keep_unmapped) out.addAlignment(rec);
				continue;
				}
			//no more record or 
			if(!buffer.isEmpty() &&
				rec!=null &&
				buffer.get(0).getReadName().equals(rec.getReadName()) &&
				buffer.get(0).getReferenceIndex().equals(rec.getReferenceIndex())
				)
				{
				buffer.add(rec);
				}
			else if(buffer.isEmpty() && rec!=null)
				{
				buffer.add(rec);
				}
			else //dump buffer
				{
				if(!buffer.isEmpty())
					{
					final int tid = buffer.get(0).getReferenceIndex();
					if(prev_tid==-1 || prev_tid!=tid)
						{
						final SAMSequenceRecord ssr= this.dict.getSequence(tid);
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
							throw new SAMException("Cigar missing in "+rec.getSAMString());
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
									final int cov = (int)readposition2coverage.incr(refPos1+x);
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
						for(final SAMRecord sr:buffer)
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
		}
	


	public static void main(final String[] args) throws IOException
		{
		new Biostar154220().instanceMainWithExit(args);
		}
	}
