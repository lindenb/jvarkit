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
package com.github.lindenb.jvarkit.tools.viewmate;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

import java.util.Set;
import java.util.function.Function;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;


/*
BEGIN_DOC

## How it works

Two modes:

  * The streaming mode is set if the input is `stdin` or if the bam file is NOT indexed. The input is scanned and any read (or mate) that overlap a region is written.

  * The other mode use the bam index. First we scan the regions, we collect the other regions and the names of the reads (requires memory), the bam is opened a second time and we collect the reads.

## Example

```
$ java -jar dist/samviewwithmate.jar -r "9:137230721-137230796"  ./src/test/resources/HG02260.transloc.chr9.14.bam | cut -f 1-9 | tail
ERR251239.10989793	83	9	137230747	60	30S70M	=	137230326	-490
ERR251239.3385449	147	9	137230754	60	1S99M	=	137230352	-500
ERR251240.17111373	99	9	137230764	60	100M	=	137231150	475
ERR251240.46859433	147	9	137230777	60	65S35M	=	137230342	-469
ERR251240.74563730	147	9	137230787	60	1S99M	=	137230407	-478
ERR251240.1291708	83	9	137230789	60	100M	=	137230411	-477
ERR251240.11887757	97	9	137230795	37	100M	14	79839451	0
ERR251239.34016218	81	14	79839349	37	100M	9	137230679	0
ERR251240.10196873	81	14	79839368	37	100M	9	137230721	0
ERR251240.11887757	145	14	79839451	37	100M	9	137230795	0
```


END_DOC
*/
@Program(name="samviewwithmate",
	description="Extract reads within given region(s), and their mates",
	keywords={"sam","bam"},
	biostars={151403,105714,368754},
	creationDate="2019-02-07",
	modificationDate="2019-02-13"
	)
public class SamViewWithMate
	extends Launcher
	{
	private static final Logger LOG = Logger.build(SamViewWithMate.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-r","--region"},description="One or more region." + IntervalParser.OPT_DESC,splitter=NoSplitter.class)
	private List<String> regionList = new ArrayList<>();
	@Parameter(names={"-b","--bed"},description="Bed file containing the region.")
	private File bedFile = null;
	@Parameter(names={"-st","--streaming"},description="Force Streaming mode even if bam is indexed. Warning: Streaming mode doesn't garantee that all mates will be fetched because a read only contains the start position of the mate of which may be out of the user's intervals, unless the MC (mate cigar) attribute is defined.")
	private boolean forceStreaming=false;
	@Parameter(names={"-u","--unmapped"},description="Also search for the unmapped mates. Not available in streaming mode.")
	private boolean look_in_unmapped = false;

	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();

	
		
	@Override
	public int doWork(final List<String> args) {
		SamReader samFileReader=null;
		SAMFileWriter sw=null;
		SAMRecordIterator  iter=null;
		try
			{
			final String input = oneFileOrNull(args);
			samFileReader = openSamReader(input);
			final SAMFileHeader header = samFileReader.getFileHeader();
			
			
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
			final Function<Locatable,QueryInterval> interval2query = RGN->{
				int tid = dict.getSequenceIndex(RGN.getContig());
				if(tid<0) throw new JvarkitException.ContigNotFoundInDictionary(RGN.getContig(), dict);
				return new QueryInterval(tid, RGN.getStart(), RGN.getEnd());
				};
			final List<QueryInterval> queryList= new ArrayList<>();
			
			final IntervalParser intervalParser=new IntervalParser(dict);
			for(final String str:regionList) {
				if(StringUtils.isBlank(str)) continue;
				final Interval interval= intervalParser.parse(str);
				if(interval==null) {
					LOG.error("Cannot parse region "+str);
					return -1;
					}
				queryList.add(interval2query.apply(interval));
				}
			if(this.bedFile!=null) {
				try (BufferedReader br=IOUtil.openFileForBufferedReading(this.bedFile)) {
					String line;
					final BedLineCodec codec=new BedLineCodec();
					final ContigNameConverter ctgConverter = ContigNameConverter.fromOneDictionary(dict);
					while((line=br.readLine())!=null) {
						final BedLine bed = codec.decode(line);
						if(bed==null) {
							LOG.warn("skipping "+line);
							continue;
							}
						final String ctg= ctgConverter.apply(bed.getContig());
						if(ctg==null) throw new JvarkitException.ContigNotFoundInDictionary(bed.getContig(), dict);
						int tid = dict.getSequenceIndex(bed.getContig());
						if(tid<0) throw new JvarkitException.ContigNotFoundInDictionary(bed.getContig(), dict);
						queryList.add(interval2query.apply(bed));
						}
					}
				}
			Collections.sort(queryList);
			
			final Function<QueryInterval,Interval> query2interval= REC->{
				return new Interval(
						dict.getSequence(REC.referenceIndex).getSequenceName(),
						REC.start,
						REC.end
						);
				};

			
			final Function<SAMRecord,Interval> mate2interval= REC->{
				if(!REC.getReadPairedFlag()) throw new IllegalStateException();
				if(REC.getMateUnmappedFlag()) throw new IllegalStateException();
			
				return new Interval(
						REC.getMateReferenceName(),
						REC.getMateAlignmentStart(),
						(SAMUtils.getMateCigar(REC)!=null?SAMUtils.getMateAlignmentEnd(REC):REC.getMateAlignmentStart())
						);
				};
			final ProgressFactory.Watcher<SAMRecord> progress=ProgressFactory.newInstance().logger(LOG).dictionary(dict).build();
			JVarkitVersion.getInstance().addMetaData(this, header);
			
			final SAMProgramRecord smr = header.createProgramRecord();
			smr.setProgramName(this.getProgramName());
			smr.setProgramVersion(this.getGitHash());
			smr.setCommandLine(this.getProgramCommandLine());
			
        	sw = this.writingBamArgs.openSAMFileWriter(this.outputFile,header, true);

			if(queryList.isEmpty())
				{
				LOG.warn("NO REGION DEFINED!!");
				// nothing to do
				}
			else if(forceStreaming || input==null || !samFileReader.hasIndex()) {
				if(!forceStreaming) {
					LOG.warning("input is stdin or bam is not indexed. Using a streaming mode");
					}
				if(look_in_unmapped)
					{
					LOG.error("cannot use --unmapped in streaming mode");
					return -1;
					}
				final IntervalTreeMap<Boolean> intervalTreeMap=new IntervalTreeMap<>();
				Arrays.stream(QueryInterval.optimizeIntervals(queryList.toArray(new QueryInterval[queryList.size()]))).
					map(query2interval).
					forEach(R->intervalTreeMap.put(R, Boolean.TRUE));
				iter=samFileReader.iterator();
				while(iter.hasNext()) {
					final SAMRecord rec = progress.apply(iter.next());

					if((!rec.getReadUnmappedFlag() && intervalTreeMap.containsOverlapping(rec)) ||
						(rec.getReadPairedFlag() && !rec.getMateUnmappedFlag() && intervalTreeMap.containsOverlapping(mate2interval.apply(rec)))) {
						rec.setAttribute(SAMTag.PG.name(), smr.getId());
						sw.addAlignment(rec);
						}
					}
				iter.close();
				iter=null;
				}
			else
				{
				boolean search_in_unmapped_flag = false;
				if(!header.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
					LOG.error("input is not sorted on coordinate.");
					return -1;
					}
				QueryInterval queryArray[]= QueryInterval.optimizeIntervals(queryList.toArray(new QueryInterval[queryList.size()]));

				final IntervalTreeMap<Boolean> intervalTreeMap=new IntervalTreeMap<>();
				Arrays.stream(queryArray).
					map(query2interval).
					forEach(R->intervalTreeMap.put(R, Boolean.TRUE));
						
				final Set<String> readNames = new HashSet<>(100_000);
				
				iter=samFileReader.query(queryArray, false);
				while(iter.hasNext()) {
					final SAMRecord rec = progress.apply(iter.next());
					if(!rec.getReadPairedFlag()) continue;
					
					if(rec.getMateUnmappedFlag()) {
						if(look_in_unmapped) {
							readNames.add(rec.getReadName());
							search_in_unmapped_flag = true;
							}
						}
					else 
						{
						final Interval rgn = mate2interval.apply(rec);
						if(intervalTreeMap.containsOverlapping(rgn)) continue;//we'll catch it anyway
						readNames.add(rec.getReadName());						
						final QueryInterval Q0= new QueryInterval(rec.getMateReferenceIndex(), rgn.getStart(),rgn.getEnd());
						queryList.add(Q0);
						}
						
					}
				iter.close();
				iter=null;
				samFileReader.close();
				
				// recompile intervals because we might have added some
				queryArray= QueryInterval.optimizeIntervals(queryList.toArray(new QueryInterval[queryList.size()]));
				// re-open			
				samFileReader = super.openSamReader(input);				
				iter=samFileReader.query(queryArray, false);
				while(iter.hasNext()) {
					final SAMRecord rec = iter.next();
					if(!(intervalTreeMap.containsOverlapping(rec) || readNames.contains(rec.getReadName()))) {
						continue;
						}
					rec.setAttribute(SAMTag.PG.name(), smr.getId());
					sw.addAlignment(rec);
					}
				iter.close();
				iter=null;
				
				//search unmapped if needed
				if(search_in_unmapped_flag) {
					iter=samFileReader.queryUnmapped();
					while(iter.hasNext()) {
						final SAMRecord rec = iter.next();
						if(!readNames.contains(rec.getReadName())) {
							continue;
							}
						rec.setAttribute(SAMTag.PG.name(), smr.getId());
						sw.addAlignment(rec);
						}
					iter.close();
					iter=null;
					}
				
				}
			progress.close();
			sw.close();
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(samFileReader);
			CloserUtil.close(sw);
			}
		}

	public static void main(final String[] args) throws Exception
		{
		new SamViewWithMate().instanceMainWithExit(args);
		}
	}
