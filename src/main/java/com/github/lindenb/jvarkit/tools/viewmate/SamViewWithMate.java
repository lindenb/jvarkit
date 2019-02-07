
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
import java.nio.file.Path;
import java.nio.file.Paths;
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

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;


/*
BEGIN_DOC

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
	biostars={151403,105714},
	creationDate="2019-02-07"
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
	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();

	
		
	@Override
	public int doWork(final List<String> args) {
		SamReader samFileReader=null;
		SAMFileWriter sw=null;
		SAMRecordIterator  iter=null;
		try
			{
			final String input = oneAndOnlyOneFile(args);
			final Path bamFile = Paths.get(input);
			IOUtil.assertFileIsReadable(bamFile);
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(bamFile);
			final Function<Locatable,QueryInterval> interval2query= RGN->{
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
			if(queryList.isEmpty()) {
				LOG.warn("NO REGION DEFINED!!");
				}
			
			QueryInterval queryArray[]= QueryInterval.optimizeIntervals(queryList.toArray(new QueryInterval[queryList.size()]));
			final IntervalTreeMap<Boolean> intervalTreeMap=new IntervalTreeMap<>();
			Arrays.stream(queryArray).
				map(A->new Interval(dict.getSequence(A.referenceIndex).getSequenceName(),A.start,A.end)).
				forEach(R->intervalTreeMap.put(R, Boolean.TRUE));
			
						
			final Set<String> readNames = new HashSet<>(100_000);
			final SamReaderFactory srf = super.createSamReaderFactory();
			samFileReader = srf.open(bamFile);
			if(!samFileReader.hasIndex()) {
				LOG.error("Sam file is not indexed "+bamFile);
				return -1;
				}
			
			
			iter=samFileReader.query(queryArray, false);
			while(iter.hasNext()) {
				final SAMRecord rec = iter.next();
				
				if(rec.getReadPairedFlag()) {
					if(!rec.getMateUnmappedFlag()) {
						final Interval rgn = new Interval(rec.getMateReferenceName(), rec.getMateAlignmentStart(), rec.getAlignmentStart());
						if(intervalTreeMap.containsOverlapping(rgn)) continue;
						readNames.add(rec.getReadName());						
						final QueryInterval Q0= new QueryInterval(rec.getMateReferenceIndex(), rgn.getStart(),rgn.getEnd());
						//LOG.debug("Adding "+rec.getReadName()+ " "+rgn);
						queryList.add(Q0);
						}
					}
				}
			iter.close();
			samFileReader.close();
			
			// recompile intervals because we might have added some
			queryArray= QueryInterval.optimizeIntervals(queryList.toArray(new QueryInterval[queryList.size()]));
			// re-open			
			samFileReader = srf.open(bamFile);
			final SAMFileHeader header = samFileReader.getFileHeader();
			JVarkitVersion.getInstance().addMetaData(this, header);
	        	sw = this.writingBamArgs.openSAMFileWriter(this.outputFile,header, true);
			iter=samFileReader.query(queryArray, false);
			while(iter.hasNext()) {
				final SAMRecord rec = iter.next();
				if(!(intervalTreeMap.containsOverlapping(rec) || readNames.contains(rec.getReadName()))) {
					continue;
					}
				sw.addAlignment(rec);
				}
			iter.close();
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
