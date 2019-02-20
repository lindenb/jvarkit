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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.pcr;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;


/**

BEGIN_DOC



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
	biostars={147136,178308},
	keywords={"sam","bam","pcr","bed"},
	modificationDate="20190220"
	)
public class PcrClipReads extends Launcher
	{
	private static final Logger LOG = Logger.build(PcrClipReads.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-B","--bed"},description="Bed file containing non-overlapping PCR fragments",required=true)
	private File bedFile = null;

	@Parameter(names={"-flag","--flag"},description="Only run on reads having sam flag (flag). -1 = all reads. (as https://github.com/lindenb/jvarkit/issues/43)")
	private int onlyFlag = -1 ;

	@Parameter(names={"-largest","--largest"},description="see if a read overlaps two bed intervals use the bed region sharing the longest sequence with a read. see https://github.com/lindenb/jvarkit/issues/44")
	private boolean chooseLargestOverlap = false;

	@Parameter(names={"-pr","--programId"},description="add a program group PG to the clipped SAM records")
	private boolean programId = false;

	
	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();
	
	private final IntervalTreeMap<Interval> bedIntervals=new IntervalTreeMap<Interval>();
	
	private Interval findInterval(final SAMRecord rec)
		{
		if(rec.getReadUnmappedFlag()) return null;
		return findInterval(rec.getContig(), rec.getAlignmentStart(), rec.getAlignmentEnd());
		}
	private Interval findInterval(final String chrom,final int start,final int end)
		{
		final Interval i= new Interval(chrom,start,end);
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
	
	
	private int run(final SamReader reader)
		{
		final SAMFileHeader header1= reader.getFileHeader();
		final SAMFileHeader header2 = header1.clone();
		Optional<SAMProgramRecord> samProgramRecord = Optional.empty();
		if(this.programId) {
			final SAMProgramRecord spr = header2.createProgramRecord();
			samProgramRecord = Optional.of(spr);
			spr.setProgramName(PcrClipReads.class.getSimpleName());
			spr.setProgramVersion(this.getGitHash());
			spr.setCommandLine(getProgramCommandLine().replace('\t', ' '));
			}
		JVarkitVersion.getInstance().addMetaData(this, header2);
		header2.setSortOrder(SortOrder.unsorted);
		SAMFileWriter sw=null;
		SAMRecordIterator iter = null;
		try
			{
			sw = this.writingBamArgs.openSAMFileWriter(this.outputFile,header2, false);
			
			final ProgressFactory.Watcher<SAMRecord> progress = ProgressFactory.newInstance().dictionary(header1).logger(LOG).build();
			iter =  reader.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec= progress.apply(iter.next());
				
				if(this.onlyFlag!=-1 &&  (rec.getFlags() & this.onlyFlag) != 0) {
					sw.addAlignment(rec);
					continue;
				}
				
				if(rec.getReadUnmappedFlag())
					{
					sw.addAlignment(rec);
					continue;
					}
				final Interval fragment = findInterval(rec);
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
				if(samProgramRecord.isPresent()) {
					readClipper.setProgramGroup(samProgramRecord.get().getId());
				}
				rec = readClipper.clip(rec, fragment);
				sw.addAlignment(rec);
				}
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
			CloserUtil.close(sw);
			}
		}
	@Override
	public int doWork(final List<String> args) {
		if(this.bedFile==null)
			{
			LOG.error("undefined bed file ");
			return -1;
			}
		BufferedReader r=null;
		SamReader samReader=null;
		try {
			final String inputName = oneFileOrNull(args);
			samReader = openSamReader(inputName);
			final SAMFileHeader header = samReader.getFileHeader();
			if(header==null)
				{
				LOG.error("No SAM header in input");
				return -1;
				}
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
			final ContigNameConverter contigNameConverter = ContigNameConverter.fromOneDictionary(dict);
			LOG.info("reading bed File "+this.bedFile);
			r= IOUtils.openFileForBufferedReading(this.bedFile);
			String line;
			final BedLineCodec codec = new BedLineCodec();
			while((line=r.readLine())!=null)
				{
				final BedLine bedLine = codec.decode(line);
				if(bedLine==null)
					{
					LOG.warn("Ignoring bed line "+line);
					continue;
					}
				Interval i = bedLine.toInterval();
				final String ctg = contigNameConverter.apply(i.getContig());
				
				if(StringUtils.isBlank(ctg)) {
					throw new JvarkitException.ContigNotFoundInDictionary(bedLine.getContig(), dict);
					}
				else if(!i.getContig().equals(ctg))
					{
					i = new Interval(ctg,i.getStart(),i.getEnd());
					}
				
				this.bedIntervals.put(i, i);
				}
			CloserUtil.close(r);r=null;
			
			return run(samReader);
			}
		catch (final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(samReader);
			this.bedIntervals.clear();;
			}
		}

	
	public static void main(final String[] args) {
		new PcrClipReads().instanceMainWithExit(args);
		}

}
