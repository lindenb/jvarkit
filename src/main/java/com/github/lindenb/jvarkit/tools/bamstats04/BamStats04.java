/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.OptionalInt;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.IntUnaryOperator;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.stats.Percentile;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;

/**
BEGIN_DOC

## input

Input is one or more indexed BAM files.

## History

* 2018-01-30: now using a jexl parser
* 2018-01-30: allow multiple values for '-cov'
* 2018-01-29: fixed bug from previous release (no data produced if no read). Added BioDas Resource.
* 2017-11-20: added new column 'partition'
* 2017-11-20: can read more than one BAM File.

## Example

```
$ java -jar dist/bamstats04.jar -B src/test/resources/toy.bed.gz src/test/resources/toy.bam 2> /dev/null | column -t 

#chrom  start  end  length  sample  mincov  maxcov  meancov  mediancov  nocoveragebp  percentcovered
ref     10     13   3       S1      3       3       3.0      3.0        0             100
ref2    1      2    1       S1      2       2       2.0      2.0        0             100
ref2    13     14   1       S1      6       6       6.0      6.0        0             100
ref2    16     17   1       S1      6       6       6.0      6.0        0             100
```

## Cited in:

  *  Han Ming Gan & al. , Genomic evidence of neo-sex chromosomes in the eastern yellow robin, GigaScience, Volume 8, Issue 9, September 2019, giz111, https://doi.org/10.1093/gigascience/giz111
  * Sigeman, H., Downing, P.A., Zhang, H. et al. The rate of W chromosome degeneration across multiple avian neo-sex chromosomes. Sci Rep 14, 16548 (2024). https://doi.org/10.1038/s41598-024-66470-7

END_DOC
 */
@Program(name="bamstats04",
	description="Coverage statistics for a BED file.",
	keywords={"sam","bam","coverage","depth","statistics","bed"},
	jvarkit_amalgamion = true,
	biostars= {309673,348251},
	creationDate = "20130513",
	modificationDate="20191003",
	menu="BAM Manipulation"
	)
public class BamStats04 extends Launcher
	{
	private static final Logger LOG = Logger.of(BamStats04.class);

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-cov","--cov"},description="add this min coverage value to ask wether the position is not covered. Use with care: any depth below this treshold will be trimmed to zero.")
	private List<Integer> minCoverages = new ArrayList<>() ;
	@Parameter(names={"-f","--filter","--jexl"},description=SamRecordJEXLFilter.FILTER_DESCRIPTION,converter=SamRecordJEXLFilter.StringConverter.class)
	private SamRecordFilter filter  = SamRecordJEXLFilter.buildDefault();
	@Parameter(names={"-B","--bed"},description=IntervalListProvider.OPT_DESC,converter=IntervalListProvider.StringConverter.class,required=true)
	private IntervalListProvider intervalListProvider = null;
	@Parameter(names={"-R","--ref"},description="[20180126]If set, a column with the GC% will be added. Also used to read CRAM. " + INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path faidxUri = null;
	@Parameter(names={"-partition","--partition"},description="[20171120]"+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition partition = SAMRecordPartition.sample;
	
	
	
	private static class IntervalStat
		{	
		private final Locatable bedLine;
		private final int counts[];
		IntervalStat(final Locatable bedLine) {
			this.bedLine = bedLine;
			this.counts=new int[bedLine.getEnd()-bedLine.getStart()+1];
			Arrays.fill(this.counts, 0);
			}
		void visit(final SAMRecord rec) {
			final Cigar cigar=rec.getCigar();
			if(cigar==null) return;
			int refpos1=rec.getAlignmentStart();
    		for(final CigarElement ce:cigar)
    			{
    			final CigarOperator op=ce.getOperator();
    			if(!op.consumesReferenceBases()) continue;
    			if(op.consumesReadBases())
    				{
    				for(int i=0;i< ce.getLength();++i)
		    			{
						if(refpos1+i>= this.bedLine.getStart() && refpos1+i<= this.bedLine.getEnd())
							{
							this.counts[refpos1+i-this.bedLine.getStart()]++;
							}
	    				}
    				}
    			refpos1+=ce.getLength();
    			if(refpos1>this.bedLine.getEnd()) break;
    			}
			}
		
		}
	
	@Override
		public int doWork(final List<String> args) {
			
			final List<Path> bamPaths = IOUtils.unrollPaths(args);
			if(bamPaths.isEmpty())  {
				LOG.error("Bam files missing");
				return -1;
				}
			
			if(this.minCoverages.isEmpty())
				{
				this.minCoverages.add(0);
				}
			
			final String NO_PARTITION="N/A";
			final List<SamReader> samReaders = new ArrayList<>(args.size());
			PrintWriter pw = null;
			ReferenceSequenceFile indexedFastaSequenceFile=null;
			GenomicSequence genomicSequence = null;
			SAMSequenceDictionary fastaDict = null;
			try
				{
				final Set<String> all_partitions = new TreeSet<>();
				SAMSequenceDictionary samDict = null;
				
				final SamReaderFactory srf = super.createSamReaderFactory();
				if(this.faidxUri!=null) srf.referenceSequence(faidxUri);
				
				for(final Path filename:bamPaths) {
					final SamReader samReader = srf.open(filename);
					if(!samReader.hasIndex()) {
						LOG.error(filename+" is not indexed");
						samReader.close();
						return -1;
						}
					final SAMFileHeader samFileheader= samReader.getFileHeader();
					if(samFileheader==null)
						{
						LOG.error("SAM file is missing a header "+filename);
						return -1;
						}
					
					final List<SAMReadGroupRecord> readGroups = samFileheader.getReadGroups();
					
					if(readGroups==null || readGroups.isEmpty())
						{
						LOG.warn("No Read group (RG) in the header of "+filename+". see https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups ");
						all_partitions.add(NO_PARTITION);
						}
					else
						{
						for(final SAMReadGroupRecord rg: readGroups)
							{
							all_partitions.add(this.partition.apply(rg,NO_PARTITION));
							}
						}
					final SAMSequenceDictionary d = SequenceDictionaryUtils.extractRequired(samFileheader);
					
					
					samReaders.add(samReader);
						
					if(samDict==null) {
						samDict=d;
						}
					else if(!SequenceUtil.areSequenceDictionariesEqual(d, samDict)) {
						LOG.error(JvarkitException.DictionariesAreNotTheSame.getMessage(d, samDict));
						return -1;
						}
					}
				
				if(samReaders.isEmpty()) {
					LOG.error("No Bam defined");
					return -1;
				}
				final ContigNameConverter samCtgConverter = ContigNameConverter.fromOneDictionary(samDict);
				
				if(this.faidxUri!=null) {
					indexedFastaSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidxUri);
					fastaDict = SequenceDictionaryUtils.extractRequired(indexedFastaSequenceFile);
					if(!SequenceUtil.areSequenceDictionariesEqual(fastaDict, samDict)) {
						LOG.error(JvarkitException.DictionariesAreNotTheSame.getMessage(fastaDict, samDict));
						return -1;
						}
					}
				pw = super.openPathOrStdoutAsPrintWriter(this.outputFile);
				pw.print(
					"#chrom\tstart\tend\tlength\t"+
					this.partition.name()+
					(indexedFastaSequenceFile==null?"":"\tgc_percent")
					);
				
				pw.print("\tmincov\tmaxcov");
				
				for(final int MIN_COVERAGE:this.minCoverages)
					{
					pw.print(
							"\tavgcov_"+MIN_COVERAGE+
							"\tmediancov_"+MIN_COVERAGE+
							"\tnocoveragebp_"+MIN_COVERAGE+
							"\tpercentcovered_"+MIN_COVERAGE
							);
					}
				pw.println();
	
			
				String line=null;
				for(final Iterator<? extends Locatable> iter2= intervalListProvider.stream().iterator();
						iter2.hasNext();
						)
					{
					final Locatable bedLine = iter2.next();
					final String ctg2 = samCtgConverter.apply(bedLine.getContig());
					if(StringUtils.isBlank(ctg2))
						{
						LOG.error("Unknown contig in "+line);
						return -1;
						}
					
					if(bedLine.getStart()>bedLine.getEnd())
						{
						LOG.info("ignoring "+bedLine);
						continue;
						}
					if(indexedFastaSequenceFile!=null && (genomicSequence==null || !genomicSequence.getChrom().equals(ctg2))) {
						if(fastaDict.getSequence(ctg2)!=null) {
							genomicSequence = new GenomicSequence(indexedFastaSequenceFile,bedLine.getContig());
							}
						else
							{
							genomicSequence = null;
							}
						}
					
					final Map<String, IntervalStat> sample2stats= new HashMap<>(all_partitions.size());
					for(final String rgId:all_partitions) {
						sample2stats.put(rgId, new IntervalStat(bedLine));
						}
					
					for(final SamReader samReader:samReaders) 
						{
						/**
						 *     start - 1-based, inclusive start of interval of interest. Zero implies start of the reference sequence.
			    		*	   end - 1-based, inclusive end of interval of interest. Zero implies end of the reference sequence. 
						 */
						final SAMRecordIterator r=samReader.queryOverlapping(
								ctg2,
								bedLine.getStart(),
								bedLine.getEnd()
								);
						while(r.hasNext())
							{
							final SAMRecord rec=r.next();
							if(rec.getReadUnmappedFlag()) continue;
							if(this.filter.filterOut(rec)) continue;
							if(!rec.getReferenceName().equals(ctg2)) continue;
							
							
							final String partition;
							final SAMReadGroupRecord group = rec.getReadGroup();
							if(group==null)
								{
								partition=NO_PARTITION;
								}
							else
								{
								final String name = this.partition.apply(group);
								partition = (StringUtil.isBlank(name)?NO_PARTITION:name);
								}
							
							IntervalStat stat= sample2stats.get(partition);
							if(stat==null) 
								{
								stat = new IntervalStat(bedLine);
								sample2stats.put(partition,stat);
								}
							stat.visit(rec);
							}
						
						r.close();
						} // end of loop over sam Readers
					
					final OptionalInt gcPercentInt = (genomicSequence==null?
						OptionalInt.empty():
						genomicSequence.getGCPercent(bedLine.getStart()-1,bedLine.getEnd()).getOptGCPercent()
						);
					
					
					for(final String partitionName : sample2stats.keySet()) {
						final IntervalStat stat = sample2stats.get(partitionName);
						Arrays.sort(stat.counts);
						
						pw.print(
								ctg2+"\t"+
								(bedLine.getStart()-1)+"\t"+
								(bedLine.getEnd())+"\t"+
								stat.counts.length+"\t"+
								partitionName
								);
						if(indexedFastaSequenceFile!=null) {
							pw.print("\t");
							if(gcPercentInt.isPresent()) pw.print(gcPercentInt.getAsInt());
							}
						pw.print(
							"\t"+
							stat.counts[0]+"\t"+
							stat.counts[stat.counts.length-1]
							);
						
						for(final int MIN_COVERAGE:this.minCoverages)
							{
							/** map depth to 0 if depth <= MIN_COVERAGE */
							final IntUnaryOperator depthAdjuster = (D)->(D<=MIN_COVERAGE?0:D);
	
							
							final int count_no_coverage=(int)Arrays.stream(stat.counts).
									filter(D-> depthAdjuster.applyAsInt(D)<=0).
									count()
									;
							
							final double mean= Percentile.average().evaluate(Arrays.stream(stat.counts).
									map(depthAdjuster)
									).getAsDouble();
							
			                final double median_depth = Percentile.median().evaluate(Arrays.stream(stat.counts).
									map(depthAdjuster)
									).getAsDouble();
			                
							
							pw.print("\t"+
									String.format("%.2f",mean)+"\t"+
									String.format("%.2f",median_depth)+"\t"+
									count_no_coverage+"\t"+
									(int)(((stat.counts.length-count_no_coverage)/(double)stat.counts.length)*100.0)
									);
							}
						pw.println();
						}
					}
				pw.flush();
				pw.close();pw=null;
				return RETURN_OK;
				}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(indexedFastaSequenceFile);
			CloserUtil.close(pw);
			CloserUtil.close(samReaders);
			}
		}
	
	public static void main(final String[] args) throws Exception
		{
		new BamStats04().instanceMainWithExit(args);
		}
}

