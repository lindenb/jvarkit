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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.FilteringSamIterator;
import htsjdk.samtools.filter.IntervalFilter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.DiscreteMedian;
import com.github.lindenb.jvarkit.math.RangeOfIntegers;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
 
/**
BEGIN_DOC

## Motivation

Because gatk is buggy: http://gatkforums.broadinstitute.org/gatk/discussion/8342/duplicate-columns-in-readlengthdistribution#latest

duplicate, secondary , supplementary, failing quality reads are discarded

## Input

input is a set of bam file or a file with suffix '.list' containing the path to the bam



## Example

```
$ java -jar dist/samreadlengthdistribution.jar -m INSERT_LENGTH src/test/resources/S*.bam 
#sample	COUNT-INSERT_LENGTH	MIN-INSERT_LENGTH	MAX-INSERT_LENGTH	AVG-INSERT_LENGTH	MEDIAN-INSERT_LENGTH	[-Inf/0[	[0/10[	[10/20[	[20/30[	[30/40[	[40/50[	[50/100[	[100/150[	[150/200[	[200/250[[250/300[	[300/350[	[350/400[	[400/450[	[450/500[	[500/1000[	[1000/Inf[
S1	999	310	696	501.31	502.00	0	0	0	0	0	0	0	0	0	003	18	115	336	527	0
S2	999	322	667	500.69	501.00	0	0	0	0	0	0	0	0	0	002	17	138	332	510	0
S3	999	322	667	500.69	501.00	0	0	0	0	0	0	0	0	0	002	17	138	332	510	0
S4	999	368	654	500.76	501.00	0	0	0	0	0	0	0	0	0	000	23	138	327	511	0
S5	999	333	641	499.78	501.00	0	0	0	0	0	0	0	0	0	004	26	126	330	513	0
```


END_DOC
 */
@Program(name="samreadlengthdistribution",
description="Sam read/insert length distribution",
keywords={"sam","bam","histogram"},
modificationDate="20200220",
creationDate="20160922"
)
public class SamReadLengthDistribution extends Launcher
	{
	private static final Logger LOG = Logger.build(SamReadLengthDistribution.class).make();
	
	private enum Method {
		SEQ_LENGTH,
		CIGAR_REF_LENGTH,
		CIGAR_PADDED_REF_LENGTH,
		INSERT_LENGTH
	}
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--groupby"},description="Group Reads by. "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition partition=SAMRecordPartition.sample;
	@Parameter(names={"-m"},description="method, how should I get the read/insert length ?")
	private Method method = Method.SEQ_LENGTH;
	@Parameter(names={"-R","--reference"},description="For reading CRAM. "+INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path faidx = null;
	@Parameter(names={"-Q","--mapq"},description="min MAPQ")
	private int mapq = 0;
	@Parameter(names={"--regions"},description="Limit analysis to this interval. "+ IntervalListProvider.OPT_DESC,splitter=NoSplitter.class,converter=IntervalListProvider.StringConverter.class)
	protected IntervalListProvider regionFiles = null;
	@Parameter(names={"-w","--windows"},description="Histogram windows. " + RangeOfIntegers.OPT_DESC,converter=RangeOfIntegers.StringConverter.class,splitter=NoSplitter.class)
	private RangeOfIntegers ranges = new RangeOfIntegers(0,10,20,30,40,50,100,150,200,250,300,350,400,450,500,1000);


	private final Map<String,DiscreteMedian<Integer>> sample2discreteMedian =new TreeMap<>();
	
	
    public SamReadLengthDistribution()
    	{
    	}		

    /** create input SAMRecord iterator */
    private CloseableIterator<SAMRecord> openSamIterator(final SamReader sr) {
    	
    	if(this.regionFiles!=null) {
    		this.regionFiles.dictionary(SequenceDictionaryUtils.extractRequired(sr.getFileHeader()));

    		if(!sr.hasIndex()) {
    			final List<Interval> L = this.regionFiles.stream().map(x->new Interval(x)).collect(Collectors.toList());
    			final SAMRecordIterator st0 = sr.iterator();
    			return new FilteringSamIterator(st0,new IntervalFilter(L, sr.getFileHeader()));
    			}
    		else
    			{
    			return sr.query(this.regionFiles.optimizedQueryIntervals(), false);
    			}
    		}
    	return sr.iterator();
    	}
    
    private void scan(final SamReader in,Path pathName) throws IOException
    	{
    	final String defName = (pathName==null?"STDIN":pathName.toString())+"#"+this.partition.name();
    	in.getFileHeader().getReadGroups().
    		stream().
    		map(RG->this.partition.apply(RG)).
    		map(S->StringUtils.isBlank(S)?defName:S).
    		filter(S->!this.sample2discreteMedian.containsKey(S)).
    		forEach(S->this.sample2discreteMedian.put(S,new DiscreteMedian<Integer>()));
    	
    	
    	final CloseableIterator<SAMRecord> iter=openSamIterator(in);
    	while(iter.hasNext())
			{
    		final SAMRecord rec = iter.next();
    		if(rec.getReadFailsVendorQualityCheckFlag()) continue;
    		if(rec.getDuplicateReadFlag()) continue;
    		if(rec.isSecondaryOrSupplementary()) continue;
    		if(!rec.getReadUnmappedFlag() && rec.getMappingQuality()<this.mapq) continue;
    		
    		final String sampleName= this.partition.getPartion(rec,defName);
    		
    		DiscreteMedian<Integer> counter = this.sample2discreteMedian.get(sampleName);
    		if(counter==null) {
    			counter =  new DiscreteMedian<>();
    			this.sample2discreteMedian.put(sampleName,counter);
    			}
    		
    		
    		final int len;
    		switch(this.method) {
    			case SEQ_LENGTH: len = rec.getReadLength(); break;
    			case CIGAR_REF_LENGTH: {
    				if(rec.getReadUnmappedFlag()) continue;
    				final Cigar c=rec.getCigar();
    				if(c==null || c.isEmpty()) continue;
    				len = c.getReferenceLength();
    				break;
    				}
    			case CIGAR_PADDED_REF_LENGTH: {
    				if(rec.getReadUnmappedFlag()) continue;
    				final Cigar c=rec.getCigar();
    				if(c==null || c.isEmpty()) continue;
    				len = c.getPaddedReferenceLength();
    				break;
    				}
    			case INSERT_LENGTH: {
    				if(rec.getReadUnmappedFlag()) continue;
    				if(!rec.getReadPairedFlag()) continue;
    				if(rec.getMateUnmappedFlag()) continue;
    				if(!rec.getContig().equals(rec.getMateReferenceName())) continue;
    				if(!rec.getFirstOfPairFlag()) continue;//ignore 2nd
    				len = Math.abs(rec.getInferredInsertSize());
    				break;
    				}
    			default: throw new IllegalStateException("unsupported method " +this.method);
    			}
    		
    		
    		
    		counter.add(len);    		
			}
    	iter.close();
    	}
    

	@Override
	public int doWork(final List<String> inputs)
		{
		SamReader r= null;
		
		try
			{
			final List<Path> args = IOUtils.unrollPaths(inputs);
			final SamReaderFactory srf = super.createSamReaderFactory();
			if(this.faidx!=null) srf.referenceSequence(this.faidx);
			
			
			if(args.isEmpty())
				{
				r = srf.open(SamInputResource.of(stdin()));
				scan(r,null);
				r.close();
				}
			else
				{
				for(final Path filename: args)
					{
					r = srf.open(SamInputResource.of(filename));
					scan(r,filename);
					r.close();
					}
				}
			try(PrintWriter out=super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				final String mtName = this.method.name();
				out.print("#"+this.partition.name()+"\tCOUNT-"+mtName+"\tMIN-"+mtName+"\tMAX-"+mtName+"\tAVG-"+mtName+"\tMEDIAN-"+mtName);
				for(final RangeOfIntegers.Range L:this.ranges.getRanges()) {
					out.print("\t");
					out.print(L.toString());
					}
				out.println();
				for(final String sample:this.sample2discreteMedian.keySet())
					{
					final DiscreteMedian<Integer> count = this.sample2discreteMedian.get(sample);
					out.print(sample);
					out.print("\t");
					out.print(count.size());
					out.print("\t");
					out.print(count.isEmpty()?".":String.valueOf((int)count.getMin().getAsDouble()));
					out.print("\t");
					out.print(count.isEmpty()?".":String.valueOf((int)count.getMax().getAsDouble()));
					out.print("\t");
					out.print(count.isEmpty()?".":String.format("%.2f",count.getAverage().getAsDouble()));
					out.print("\t");
					out.print(count.isEmpty()?".":String.format("%.2f",count.getMedian().getAsDouble()));
					for(final RangeOfIntegers.Range L:this.ranges.getRanges()) {
						out.print("\t");
						out.print(count.getAsDoubleStream().filter(I->L.contains((int)I)).count());
						}
					out.println();
					}
				out.flush();
				}
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	
	public static void main(final String[] args) {
		new SamReadLengthDistribution().instanceMainWithExit(args);

	}

}
