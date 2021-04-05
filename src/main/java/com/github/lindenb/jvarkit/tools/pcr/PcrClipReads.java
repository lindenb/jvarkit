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
package com.github.lindenb.jvarkit.tools.pcr;

import java.util.ArrayList;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassBamLauncher;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;


/**

BEGIN_DOC


## See also

* (2021) `samtools ampliconclip` – clip reads using a BED file  http://www.htslib.org/doc/samtools-ampliconclip.html


## Motivation


 Soft clip BAM files based on PCR target regions https://www.biostars.org/p/147136/


 *  mapping quality is set to zero if a read on mapped strand - overlap the 5' side of the PCR fragment
 *  mapping quality is set to zero if a read on mapped strand + overlap the 3' side of the PCR fragment
 *  mapping quality is set to zero if no PCR fragment is found


after processing the BAM file should be sorted, processed with samtools calmd and picard fixmate


### Example


```
echo  "seq2\t1100\t1200" > test.bed
java -jar dist/pcrclipreads.jar -B test.bed  samtools-0.1.19/examples/ex1.bam  |\
	samtools  view -q 1 -F 4 -Sbu  -  |\
	samtools  sort -o clipped.bam -  && samtools index clipped.bam

samtools tview -p seq2:1100  clipped.bam  samtools-0.1.19/examples/ex1.fa

```


### output


![img](http://i.imgur.com/bjDEnMW.jpg)



```
    1091      1101      1111      1121      1131      1141      1151      1161      1171      1181      1191
AAACAAAGGAGGTCATCATACAATGATAAAAAGATCAATTCAGCAAGAAGATATAACCATCCTACTAAATACATATGCACCTAACACAAGACTACCCAGATTCATAAAACAAATNNNNN
              ...................................                               ..................................
              ,,,                                                               ..................................
              ,,,,,                                                              .................................
              ,,,,,,,,,,,                                                        .............................N...
              ,,,,,,,,                                                             ...............................
              ,,g,,,,,,,,,,,,,                                                        ............................
              ,,,,,,,,,,,,,,,,,,,,                                                    ............................
              ,,,,,,,,,,,,,,,,,,,                                                       ..........................
              ,,,,,,,,,,,,,,,,,,,,,,                                                    ..........................
              ,,,,,,,,,,,,,,,,,,,,,,,                                                       ......................
              ,,,,,,,,,,,,,,,,,,,,,,,,,,                                                        ..................
              ,,,,,,,,,,,,,,,,,,,,,,,,,,                                                        ..................
              ,,,,,,,,,,,,,,,,,,,,,,,,,,,,                                                       .................
              ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,                                                       ................
              ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,                                                        ...............
              ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,                                                         ............
              ,,,,,,,,,,,,a,,,,,,,,,,,,,,,,,,,                                                             .......
              ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,                                                            ......
              ,,a,,,a,,,,,,,,,,,,,,,,,,,,,,,,,,,                                                              ....
              ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,                                                             ....
              ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,                                                                .
                                                                                                                 .

```

## Cited in

 * BAMClipper: removing primers from alignments to minimize false-negative mutations in amplicon next-generation sequencing [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5431517/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5431517/)


## History

 * 20170630 : rewritten after [https://github.com/lindenb/jvarkit/issues/81](https://github.com/lindenb/jvarkit/issues/81)



END_DOC
*/

@Program(name="pcrclipreads",
	description="Soft clip bam files based on PCR target regions",
	biostars={147136,178308,498088},
	keywords={"sam","bam","pcr","bed"},
	modificationDate="20210322",
	creationDate="20150618"
    )
public class PcrClipReads extends OnePassBamLauncher {
	private static final Logger LOG = Logger.build(PcrClipReads.class).make();


	@Parameter(names={"-B","--bed","--pcr"},description="Regions containing non-overlapping PCR fragments. "+IntervalListProvider.OPT_DESC,converter=IntervalListProvider.StringConverter.class,required=true)
	private IntervalListProvider intervalListProvider = IntervalListProvider.unspecified();

	@Parameter(names={"-flag","--flag"},description="Only run on reads having sam flag 'x' (flag). -1 = all reads. (as https://github.com/lindenb/jvarkit/issues/43)")
	private int onlyFlag = -1 ;

	@Parameter(names={"-largest","--largest"},description="check if a read overlaps two bed intervals use the bed region sharing the longest sequence with a read. see https://github.com/lindenb/jvarkit/issues/44")
	private boolean chooseLargestOverlap = false;

	@Parameter(names={"-pr","--programId"},description="add a program group PG to the clipped SAM records")
	private boolean programId = false;
	
	private final IntervalTreeMap<Interval> bedIntervals=new IntervalTreeMap<Interval>();
	private SAMProgramRecord samProgramRecord=null;
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	
	private Interval findInterval(final SAMRecord rec)
		{
		if(rec.getReadUnmappedFlag()) return null;
		return findInterval(new Interval(rec));
		}
	private Interval findInterval(final Interval i)
		{
		final List<Interval> L= new ArrayList<>(this.bedIntervals.getOverlapping(i));
		if(L.isEmpty()) return null;
		if(L.size()==1) return L.get(0);
		if(this.chooseLargestOverlap)
			{
			Interval best = null;
			Integer dist = null;
			for(final Interval j: L) {
				if(!i.intersects(j)) continue;//????
				final int commonDist = i.intersect(j).length();
				if(dist==null || dist < commonDist) {
					best = j;
					dist = commonDist;
					}	
				}
			return best;
			}
		else
			{
			throw new IllegalStateException(
					"Option chooseLargestOverlap not used. Overlapping PCR intervals samRecord "+i+": "+L);
			}
		}
	
	@Override
	protected SAMFileHeader createOutputHeader(SAMFileHeader headerIn) {
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(headerIn);
		this.intervalListProvider.
			dictionary(dict).
			stream().
			map(R->new Interval(R)).
			forEach(R->
			{
			this.bedIntervals.put(R,R);
			});
		
		final SAMFileHeader header2 = super.createOutputHeader(headerIn);
	

		if(this.programId) {
			this.samProgramRecord = super.createProgramRecord(header2);
			}
		header2.setSortOrder(SortOrder.unsorted);
		return header2;
		}
	
	@Override
	protected void scanIterator(final SAMFileHeader headerIn, CloseableIterator<SAMRecord> iter, SAMFileWriter sw) {
		while(iter.hasNext())
			{
			SAMRecord rec= iter.next();
			
			if(this.onlyFlag!=-1 &&  (rec.getFlags() & this.onlyFlag) != 0) {
				sw.addAlignment(rec);
				continue;
				}
			
			if(rec.getReadUnmappedFlag())
				{
				sw.addAlignment(rec);
				continue;
				}
			final Interval fragment = this.findInterval(rec);
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
			//    REC     >>>>>>>
			//    FRAG       xxxxxxxxx
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
			final ReadClipper readClipper = new ReadClipper();
			if(this.samProgramRecord!=null) {
				readClipper.setProgramGroup(this.samProgramRecord.getId());
				}
			rec = readClipper.clip(rec, fragment);
			sw.addAlignment(rec);
			}// end while
		}	
		
	
	
	
	public static void main(final String[] args) {
		new PcrClipReads().instanceMainWithExit(args);
		}

}
