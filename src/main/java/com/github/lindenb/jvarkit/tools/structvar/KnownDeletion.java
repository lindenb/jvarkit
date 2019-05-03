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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamRecordFilterFactory;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;


/**
BEGIN_DOC

## Motivation

search known deletion
## Example

```
find DIR -type f -name "*.bam" > bam.list

```

END_DOC

 */
@Program(name="knowndel",
	description="Find Split reads in a region to validate a known CNV",
	keywords= {"cnv","bam","sam","vcf","sv","deletion"},
	creationDate="20190523",
	modificationDate="20190523",
	generate_doc=false
	)
public class KnownDeletion extends Launcher
	{
	private static final Logger LOG = Logger.build(KnownDeletion.class).make();
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"-r","--rgn","--region"},description="region chr:start-end to be validated.",required=true)
	private String regionStr="";
	@Parameter(names={"-x","--extend"},description="search in the boundaries of the CNV around +/- 'x' bases. "+ DistanceParser.OPT_DESCRIPTION, converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int extendDistance = 100;	
	@Parameter(names={"-filter","--filter"},description=SamRecordFilterFactory.FILTER_DESCRIPTION,converter=SamRecordFilterFactory.class)
	private SamRecordFilter samRecordFilter = SamRecordFilterFactory.getDefault();
	@Parameter(names={"-B","--bam"},description="Optional Bam to sasve the matching reads")
	private Path bamOut=null;
	
	@Override
	public int doWork(final List<String> args) {		
		if(this.extendDistance<=0)
			{
			LOG.error("bad extend factor "+this.extendDistance);
			return -1;
			}
		
		final Function<SAMRecord, Integer> mateEnd = REC->SAMUtils.getMateCigar(REC)!=null?
				SAMUtils.getMateAlignmentEnd(REC):
				REC.getMateAlignmentStart()
				;
		
		PrintWriter pw = null;
		SAMFileWriter sfw = null;
		try
			{
			pw = super.openPathOrStdoutAsPrintWriter(this.outputFile);
			final SamReaderFactory srf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
			
			final List<String> filenames=IOUtils.unrollStrings2018(args);
			if(this.bamOut!=null && !filenames.isEmpty()) {
				SAMSequenceDictionary theDict = null;
				final Set<String> samples = new TreeSet<>();
				
				for(final String filename: filenames) {
					final SamInputResource sri = SamInputResource.of(filename);
					try(final SamReader samReader = srf.open(sri) ) {
					final SAMFileHeader header  = samReader.getFileHeader();
					final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
					if(theDict==null) {
						theDict=dict;
						}
					else if(!SequenceUtil.areSequenceDictionariesEqual(theDict, dict)) {
						throw new JvarkitException.DictionariesAreNotTheSame(theDict, dict);
						}
					final String sampleName = header.getReadGroups().stream().
								map(R->R.getSample()).
								filter(S->!StringUtil.isBlank(S)).
								findFirst().
								orElse(filename);
					samples.add(sampleName);
					}
				}
				final SAMFileHeader header=new SAMFileHeader(theDict);
				header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
				JVarkitVersion.getInstance().addMetaData(this, header);
				for(final String sn:samples) {
					final SAMReadGroupRecord srg = new SAMReadGroupRecord(sn);
					srg.setSample(sn);
					header.addReadGroup(srg);
					}
				final SAMFileWriterFactory swf = 	new SAMFileWriterFactory();
				sfw = swf.makeSAMOrBAMWriter(header, true, this.bamOut);
				}
			
			for(final String filename: IOUtils.unrollStrings2018(args)) {
				final SamInputResource sri = SamInputResource.of(filename);
				try(final SamReader samReader = srf.open(sri) ) {
				final SAMFileHeader header  = samReader.getFileHeader();
				final SAMSequenceDictionary dict = header.getSequenceDictionary();
				final String sampleName = header.getReadGroups().stream().
							map(R->R.getSample()).
							filter(S->!StringUtil.isBlank(S)).
							findFirst().
							orElse(filename);
				final IntervalParser parser = new IntervalParser(dict);
				final Interval interval =  parser.parse(this.regionStr);
				if(interval==null) {
					LOG.error("Cannot parser "+this.regionStr+" for "+filename);
					return -1;
					}
				final SAMSequenceRecord ssr = Objects.requireNonNull(dict.getSequence(interval.getContig()));
				final int tid= ssr.getSequenceIndex();
				final QueryInterval qi1 = new QueryInterval(
						tid,
						Math.max(0,interval.getStart() - this.extendDistance),
						Math.min(interval.getStart() + this.extendDistance,ssr.getSequenceLength())
						);
				final QueryInterval qi2 = new QueryInterval(
						tid,
						Math.max(0,interval.getEnd() - this.extendDistance),
						Math.min(interval.getEnd() + this.extendDistance,ssr.getSequenceLength())
						);
				if(CoordMath.overlaps(qi1.start, qi1.end, qi2.start, qi2.end)) {
					LOG.error("after extending the boundaries with "+this.extendDistance+". They both overlap...");
					return -1;
					}
				long count_disc = 0L;
				long count_split = 0L;
				long count_del = 0L;
				
				final QueryInterval intervals[] = QueryInterval.optimizeIntervals(new QueryInterval[]{qi1,qi2});
				try(final CloseableIterator<SAMRecord> iter = samReader.query(intervals, false)) {
					while(iter.hasNext()) {
						final SAMRecord rec  = iter.next();
						if(rec.getReadUnmappedFlag()) continue;
						if(this.samRecordFilter.filterOut(rec)) continue;
						
						if(rec.getStart() <= qi1.end && rec.getEnd()>=qi2.start) {
							count_del++;
							if(sfw!=null) {
								rec.setAttribute("RG",sampleName);
								sfw.addAlignment(rec);
								}
							continue;
							}
						
						
						if(rec.getReadPairedFlag() &&
							!rec.getMateUnmappedFlag() &&
							rec.getMateReferenceIndex().equals(tid) &&
							(
								( CoordMath.overlaps(rec.getStart(),rec.getEnd(), qi1.start,qi1.end) &&
								  CoordMath.overlaps(rec.getMateAlignmentStart(),mateEnd.apply(rec),qi2.start,qi2.end)) ||
								
								( CoordMath.overlaps(rec.getStart(),rec.getEnd(), qi2.start,qi2.end) &&
								  CoordMath.overlaps(rec.getMateAlignmentStart(),mateEnd.apply(rec),qi1.start,qi1.end)) 
							))
							{
							count_disc++;
							if(sfw!=null) {
								rec.setAttribute("RG",sampleName);
								sfw.addAlignment(rec);
								}
							continue;
							}
						
						for(final SAMRecord rec2:SAMUtils.getOtherCanonicalAlignments(rec))
							{
							if(	rec2.getReferenceIndex().equals(tid) &&
								(
								( CoordMath.overlaps(rec.getStart(),rec.getEnd(), qi1.start,qi1.end) && 
								  CoordMath.overlaps(rec2.getStart(),rec2.getEnd(), qi2.start,qi2.end)) ||
								
								( CoordMath.overlaps(rec.getStart(),rec.getEnd(), qi2.start,qi2.end) &&
								  CoordMath.overlaps(rec2.getStart(),rec2.getEnd(), qi1.start,qi1.end))
								))
								{
								count_split++;
								if(sfw!=null) {
									rec.setAttribute("RG",sampleName);
									sfw.addAlignment(rec);
									}
								break;
								}
							}
						}
					}
				pw.println(filename+"\t"+sampleName+"\t"+this.regionStr+"\t"+count_disc+"\t"+count_split+"\t"+count_del+"\t"+(count_del+count_disc+count_split));
				}
			}
			pw.flush();
			if(sfw!=null) sfw.close(); sfw=null;
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(sfw);
			}
		}
	
	
	public static void main(final String[] args) {
		new KnownDeletion().instanceMainWithExit(args);
		}
	}
