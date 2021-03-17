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
package com.github.lindenb.jvarkit.tools.bamstats04;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.DiscreteMedian;
import com.github.lindenb.jvarkit.samtools.CoverageFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;

/**
BEGIN_DOC

## Input 

input is one or more indexed bam file.

One file with  the suffix '.list' is interpreted as a text file with one path per line.

If there is no argument, stdin is interpreted as a list of path to the bam like in `find . -name "*.bam"`


## Cited In:

  * "Custom hereditary breast cancer gene panel selectively amplifies target genes for reliable variant calling" . BioRxiv https://doi.org/10.1101/322180
  * "Future MicrobiologyVol. 15, No. 15  Comparative genomics of Sporothrix species and identification of putative pathogenic-gene determinants.  Published Online:12 Nov 2020https://doi.org/10.2217/fmb-2019-0302"


## Example

```
$ head genes.bed
1	179655424	179655582	ZORG
1	179656788	179656934	ZORG

$ java -jar  dist/bamstats05.jar -B genes.bed --mincoverage 10 in.bam > out.txt

$ head out.txt
#chrom	start	end	gene	sample	length	mincov	maxcov	avg	nocoverage.bp	percentcovered
1	179655424	179656934	ZORG	SAMPLE1	304	27	405	216.80921052631578	0	100
```


END_DOC


*/
@Program(name="bamstats05",
description="Coverage statistics for a BED file, group by gene",
keywords={"bam","coverage","statistics","bed"},
biostars={324639,194393,35083},
modificationDate="20210317",
creationDate="20151012"
)
public class BamStats05 extends Launcher
	{
	private static final Logger LOG = Logger.build(BamStats05.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-m","--mincoverage"},description="Coverage treshold. Any depth under this value will be considered as 'not-covered'.  Default: 0")
	private List<Integer> min_coverages = new ArrayList<>() ;
	@Parameter(names={"-merge","--merge"},description="[20181122] Merge overlapping intervals for the same gene.")
	private boolean mergeOverlapping =  false;
	@Parameter(names={"-B","--bed"},description="bed file (columns: chrom(tab)start(tab)end(tab)GENE)",required=true)
	private Path BEDILE = null;
	@Parameter(names={"--groupby"},description="Group Reads by. "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition groupBy=SAMRecordPartition.sample;
	@Parameter(names={"-f","--filter","--jexl"},description=SamRecordJEXLFilter.FILTER_DESCRIPTION,converter=SamRecordJEXLFilter.StringConverter.class)
	private SamRecordFilter filter  = SamRecordJEXLFilter.buildDefault();
	@Parameter(names={"-R","--reference"},description=CRAM_INDEXED_REFENCE)
	private Path faidx = null;
	@Parameter(names={"--mapq"},description="Min mapping quality")
	private int mapq = 1;

	
	private Map<String, List<SimpleInterval>> readBedFile(final Path bedFile) throws IOException
    	{
    	final Map<String, List<SimpleInterval>> gene2interval=new TreeMap<String, List<SimpleInterval>>();
    	
    	BufferedReader bedIn=null;
    	try
    		{
    		bedIn=IOUtils.openPathForBufferedReading(bedFile);
    		final BedLineCodec codec = new BedLineCodec();
    		String line=null;
			while((line=bedIn.readLine())!=null)
				{
				if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
				final BedLine bedLine = codec.decode(line);
				if(bedLine==null) continue;
				if(bedLine.getColumnCount()<4)
					{
					throw new IOException("bad bed line  (expected 4 columns CHROM/START/END/GENE) in "+line+" "+bedFile);
					}
				final String chrom = bedLine.getContig();
				final int chromStart1= bedLine.getStart();
				final int chromEnd1= bedLine.getEnd();
				final String gene = bedLine.get(3);
				if(StringUtils.isBlank(gene))  throw new IOException("bad bed gene in "+line+" "+bedFile);
				 List<SimpleInterval> intervals = gene2interval.get(gene);
				 if(intervals==null)
				 	{
					 intervals=new ArrayList<>();
					 gene2interval.put(gene,intervals);
				 	} 
				 else if(!intervals.get(0).getContig().equals(chrom))
				 	{
					throw new IOException("more than one chromosome for gene:"+gene+" "+intervals.get(0)+" and "+bedLine.toInterval());
				 	}
				else
				 	{
					if(!this.mergeOverlapping) {
						intervals.stream().
							filter(R->R.overlaps(bedLine)).
							findFirst().
							ifPresent(R->{
								throw new IllegalArgumentException("overlapping region: "+bedLine+" and "+R+
										"\nUse option --merge to merge overlapping regions.");
								});
						}
				 	}
				intervals.add(new SimpleInterval(chrom, chromStart1, chromEnd1));
				if(this.mergeOverlapping)
					{
					intervals.sort((A,B)->Integer.compare(A.getStart(),B.getStart()));
					int x = 0;
					while(x+1<intervals.size())
						{
						final SimpleInterval i1 = intervals.get(x+0);
						final SimpleInterval i2 = intervals.get(x+1);
						if(i1.overlaps(i2))
							{
							intervals.remove(x+1);
							intervals.set(x+0,i1.merge(i2));
							}
						else
							{
							x++;
							}
						}
					}
				}
    		bedIn.close();
    		return gene2interval;
    		}
    	finally
    		{
    		CloserUtil.close(bedIn);
    		}
    	}
	
	protected  int doWork(
			final PrintWriter pw,
			final Map<String, List<SimpleInterval>> gene2interval,
			final String filename,
			final SamReader IN) throws Exception
		{
		try
			{
			LOG.info("Scanning "+filename);
			final SAMFileHeader header = IN.getFileHeader();
			final List<SAMReadGroupRecord> rgs = header.getReadGroups();
			if(rgs==null || rgs.isEmpty())
				throw new IOException("No read groups in "+filename);
			final Set<String> groupNames = this.groupBy.getPartitions(rgs);
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
			final ContigNameConverter contigNameConverter = ContigNameConverter.fromOneDictionary(dict);
			final CoverageFactory coverageFactory = new CoverageFactory().
					setMappingQuality(this.mapq).
					setPartition(this.groupBy).
					setRecordFilter(R->!filter.filterOut(R));
			
			for(final String partition : groupNames)
				{
				if(StringUtils.isBlank(partition)) throw new IOException("Empty read group: "+groupBy.name()+" for "+filename);
				for(final String gene: gene2interval.keySet())
					{
					
					final List<Integer> counts = new ArrayList<>();
					final List<SimpleInterval> intervals = gene2interval.get(gene);
					final String newContig = contigNameConverter.apply(intervals.get(0).getContig());
					if(StringUtil.isBlank(newContig)) {
						throw new JvarkitException.ContigNotFoundInDictionary(intervals.get(0).getContig(), dict);
						}
					
					final CoverageFactory.SimpleCoverage coverage = coverageFactory.getSimpleCoverage(
							IN,
							intervals.stream().map(R->R.renameContig(newContig)).collect(Collectors.toList()),
							partition
							);
					
					for(final SimpleInterval interval:intervals)
						{
						for(int i=interval.getStart();i<=interval.getEnd() && i <= coverage.getEnd();i++) {
							final int d = coverage.get(i-coverage.getStart());
							counts.add(d);
							}
						}
					
					Collections.sort(counts);
						
						
					pw.print(
							intervals.get(0).getContig()+"\t"+
							(coverage.getStart()-1)+"\t"+
							coverage.getEnd() +"\t"+gene+"\t"+partition+"\t"+
							intervals.size()+"\t"+
							counts.size()+"\t"+
							counts.get(0)+"\t"+
							counts.get(counts.size()-1)
							);
					
					for(final int mc:this.min_coverages)
						{
						final DiscreteMedian<Integer> discreteMedian = new DiscreteMedian<>();
						int count_no_coverage=0;
						for(int cov:counts)
							{
							if(cov<=mc) ++count_no_coverage;
							discreteMedian.add(cov);
							}
						
						final OptionalDouble average = discreteMedian.getAverage();
						final OptionalDouble median = discreteMedian.getMedian();
						
						
						pw.print("\t"+
								(average.isPresent()?String.format("%.2f",average.orElse(0.0)):".")+"\t"+
								(median.isPresent()?String.format("%.2f",median.orElse(-1.0)):".")+"\t"+
								count_no_coverage+"\t"+
								(int)(((counts.size()-count_no_coverage)/(double)counts.size())*100.0)
								);
						}
					
					pw.println();
					}//end gene
				}//end sample
		return 0;
		}
	catch(final Throwable err)
		{
		LOG.error(err);
		return -1;
		}
	finally
		{
		CloserUtil.close(IN);
		}
	}
	
	@Override
	public int doWork(final List<String> args) {
		if(BEDILE==null)
			{
			LOG.error( "missing bed file");
			return -1;
			}
		if(this.min_coverages.isEmpty()) min_coverages.add(0);

		try
			{
			final Map<String, List<SimpleInterval>> gene2interval = readBedFile(BEDILE);
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				//print header
				pw.print(
						"#chrom\t"+
						"gene.start.0"+"\t"+"gene.end.0"+"\t"+"gene.Name"+"\t"+groupBy.name()+"\t"+
						"count.intervals\t"+
						"length"+"\t"+
						"min.cov"+"\t"+
						"max.cov");
				
				for(final int mc:this.min_coverages)
					{
					pw.print("\t"+
							"mean.GT_"+mc+"\t"+
							"median.GT_"+mc+"\t"+
							"no_coverage.GT_"+mc+"\t"+
							"percent_covered.GT_"+mc
							);
					}		
				pw.println();
				
				final SamReaderFactory srf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
				if(this.faidx!=null) srf.referenceSequence(this.faidx);
				
				final List<Path> files = IOUtils.unrollPaths(args);
				
				
				for(final Path f:files)
					{
					try(SamReader in = srf.open(f)) {
						int tl = doWork(pw,gene2interval,f.toString(),in);
						if(tl!=0) return tl;
						}
					}
				pw.flush();
				}
			return 0;
			}
		catch (final Throwable e) {
			LOG.error(e);
			return -1;
			}
		finally
			{
			}
		}

	
	public static void main(final String[] args) throws IOException
		{
		new BamStats05().instanceMainWithExit(args);
		}

	}

