
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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.bamstats01;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;

import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;

/**

BEGIN_DOC


### Example


```
$ java -jar dist/samstats01.jar --bed "RF03:1-3000"  src/test/resources/S*.bam 
#Filename	Sample	UNMAPPED	NOT_PRIMARY	FAIL_VENDOR_QUALITY	OFF_TARGET	DUPLICATE	FAIL_MAPPING_QUALITY	OK_CALLING
src/test/resources/S1.bam	S1	0	0	0	1718	0	0280
src/test/resources/S2.bam	S2	0	0	0	1718	0	0280
src/test/resources/S3.bam	S3	0	0	0	1718	0	0280
src/test/resources/S4.bam	S4	0	0	0	1718	0	0280
src/test/resources/S5.bam	S5	0	0	0	1718	0	0280
```


END_DOC
*/


@Program(name="samstats01",
	description="Statistics about the reads in a BAM.",
	keywords= {"sam","bam"},
	modificationDate="20191004",
	creationDate="20130705"
	)
public class BamStats01
	extends Launcher
	{

	private static final Logger LOG = Logger.build(BamStats01.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;

	@Parameter(names={"-B","--bed","--regions"},description=IntervalListProvider.OPT_DESC,converter=IntervalListProvider.StringConverter.class)
	private IntervalListProvider intervalListProvider = null;
	@Parameter(names={"-q","--qual"},description="min mapping quality")
	private int minMappingQuality = 30 ;
	@Parameter(names={"--groupby"},description="Group Reads by. "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition groupBy=SAMRecordPartition.sample;
	@Parameter(names={"-R","--reference"},description=CRAM_INDEXED_REFENCE)
	private Path faidx=null;
    

	private enum Category
		{
		UNMAPPED,
		NOT_PRIMARY,
		FAIL_VENDOR_QUALITY,
		OFF_TARGET,
		DUPLICATE,
		FAIL_MAPPING_QUALITY,
		OK_CALLING
		};
	
		
	private void run(final PrintWriter out,final String filename,final SamReader samFileReader) throws IOException
		{	
		final Map<String,Counter<Category>> sample2hist=new HashMap<>();
		final SAMFileHeader header = samFileReader.getFileHeader();
		final SAMSequenceDictionary currDict=SequenceDictionaryUtils.extractRequired(header);
		final IntervalTreeMap<Locatable> intervalTreeMap;
		
		header.getReadGroups().stream().
			map(RG->groupBy.apply(RG)).
			filter(S->!StringUtils.isBlank(S)).
			map(S->sample2hist.put(S,new Counter<Category>()));
		
		
		if(this.intervalListProvider!=null) {
				intervalTreeMap = this.intervalListProvider.
						dictionary(currDict).
						toIntervalTreeMap()
						;
				}
		else {
			intervalTreeMap = null;
			}
		try(final SAMRecordIterator iter=samFileReader.iterator()) {
		while(iter.hasNext())
			{
			final SAMRecord rec= iter.next();
			
			final String sampleName = groupBy.getPartion(rec,"undefined");
			
			Counter<Category> hist=sample2hist.get(sampleName);
			if(hist==null)
				{
				hist=new Counter<>();
				sample2hist.put(sampleName, hist);
				}
			
			if(rec.getReadUnmappedFlag()) {
				hist.incr(Category.UNMAPPED);
				continue;
				}
			
			if(rec.isSecondaryOrSupplementary()) {
				hist.incr(Category.NOT_PRIMARY);
				continue;
			}
			
			if(rec.getReadFailsVendorQualityCheckFlag()) {
				hist.incr(Category.FAIL_VENDOR_QUALITY);
				continue;
			}
			
			
			if(intervalTreeMap!=null &&
				!intervalTreeMap.containsOverlapping(rec))
				{
				hist.incr(Category.OFF_TARGET);
				continue;
				}		
			
			if(rec.getMappingQuality() < this.minMappingQuality) {
				hist.incr(Category.FAIL_MAPPING_QUALITY);
				continue;
			}
			hist.incr(Category.OK_CALLING);
				
			}
		
		
		for(final String sampleName: sample2hist.keySet()) {
			final Counter<Category> hist=sample2hist.get(sampleName);
			out.print(filename+"\t"+sampleName);
			
			for(final Category cat1: Category.values())//je je suis libertineuuh, je suis une cat1
				{
				out.print("\t");
				out.print(hist.count(cat1));
				}
			}
		out.println();
		}
	out.flush();
	}
	
	@Override
	public int doWork(final List<String> inputs) {
		final List<Path> args=IOUtils.unrollPaths(inputs);
		try {
			try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				out.println("#Filename\tSample\t"+ Arrays.stream(Category.values()).map(E->E.name()).collect(Collectors.joining("\t")));	
				final SamReaderFactory srf= super.createSamReaderFactory();
				if(this.faidx!=null) {
					srf.referenceSequence(this.faidx);
					}
				
				if(args.isEmpty())
					{
					try(final SamReader r= srf.open(SamInputResource.of(stdin()))) {
						run(out,"stdin",r);
						}
					}
				else
					{
					for(final Path filename:args)
						{
						try(final SamReader sfr=srf.open(filename))
							{
							run(out,filename.toString(),sfr);
							}
						}
					}
				out.flush();			
				}
			return 0;
		} catch (final Throwable e) {
			LOG.error(e);
			return -1;
		} finally
			{
			}
		}
	
	
	public static void main(final String[] args)
		{
		new BamStats01().instanceMainWithExit(args);
		}
	}
