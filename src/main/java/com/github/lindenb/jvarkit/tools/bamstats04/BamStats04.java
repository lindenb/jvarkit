/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.bamstats04;

import java.io.BufferedReader;
import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamFilterParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

/**
BEGIN_DOC

## Example

```
$ java -jar dist/bamstats04.jar \
	-B data.bed \
	f.bam
#chrom	start	end	length	mincov	maxcov	mean	nocoveragebp	percentcovered
1	429665	429785	120	42	105	72.36666666666666	0	100
1	430108	430144	36	9	9	9.0	0	100
1	439811	439904	93	0	36	3.6451612903225805	21	77
1	550198	550246	48	1325	1358	1344.4583333333333	0	100
1	629855	629906	51	223	520	420.70588235294116	0	100
1	689960	690029	69	926	1413	1248.9420289855072	0	100
1	690852	690972	120	126	193	171.24166666666667	0	100
1	787283	787406	123	212	489	333.9756097560976	0	100
1	789740	789877	137	245	688	528.6715328467153	0	1
```

END_DOC
 */
@Program(name="bamstats04",
	description="Coverage statistics for a BED file.",
	keywords={"bam","coverage","statistics","bed"}
	)
public class BamStats04 extends Launcher
	{
	private static final Logger LOG = Logger.build(BamStats04.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-cov","--cov"},description="min coverage to say the position is not covered")
	private int MIN_COVERAGE = 0 ;

	@Parameter(names={"-f","--filter"},description=SamFilterParser.FILTER_DESCRIPTION,converter=SamFilterParser.StringConverter.class)
	private SamRecordFilter filter  = SamFilterParser.buildDefault();

	@Parameter(names={"-B","--bed"},description="Bed File. Required",required=true)
	private File bedFile = null;

	@Parameter(names={"-R","--ref"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION+" If set, a column with the GC% will be added")
	private File faidxFile = null;

	@Override
		public int doWork(final List<String> args) {
			if(this.bedFile==null || !this.bedFile.exists()) {
				LOG.error("undefined option -B");
				return -1;
			}
			BufferedReader bedIn=null;
			SamReader samReader = null;
			PrintWriter pw = null;
			IndexedFastaSequenceFile indexedFastaSequenceFile=null;
			GenomicSequence genomicSequence=null;
			try
				{
				
				final BedLineCodec codec= new BedLineCodec();
				
				bedIn=IOUtils.openFileForBufferedReading(this.bedFile);
				samReader = super.openSamReader(oneFileOrNull(args));
				final SAMSequenceDictionary dict = samReader.getFileHeader().getSequenceDictionary();
					if(dict==null) {
					LOG.error("SAM sequence dictionary missing ");
					return -1;
					}
				
				
				if(this.faidxFile!=null) {
					indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.faidxFile);
				}
				pw = super.openFileOrStdoutAsPrintWriter(this.outputFile);
				pw.println("#chrom\tstart\tend\tlength\t"+
					(indexedFastaSequenceFile==null?"":"gc_percent\t")+
					"mincov\tmaxcov\tmeancov\tmediancov\tnocoveragebp\tpercentcovered");
	
			
				String line=null;
				while((line=bedIn.readLine())!=null)
					{
					if(line.isEmpty() || line.startsWith("#")) continue;
					final BedLine bedLine = codec.decode(line);
					if(bedLine==null) continue;
					if(dict.getSequence(bedLine.getContig())==null)
						{
						LOG.error("Unknown contig in "+line);
						return -1;
						}
					
					if(indexedFastaSequenceFile!=null && (genomicSequence==null || !genomicSequence.getChrom().equals(bedLine.getContig()))) {
						genomicSequence = new GenomicSequence(indexedFastaSequenceFile, bedLine.getContig());
						}
					
					/* picard javadoc:  - Sequence name - Start position (1-based) - End position (1-based, end inclusive)  */
					
					
					final int counts[]=new int[bedLine.getEnd()-bedLine.getStart()+1];
					if(counts.length==0) continue;
					Arrays.fill(counts, 0);
					
					
					/**
					 *     start - 1-based, inclusive start of interval of interest. Zero implies start of the reference sequence.
		    		*	   end - 1-based, inclusive end of interval of interest. Zero implies end of the reference sequence. 
					 */
					final SAMRecordIterator r=samReader.queryOverlapping(
							bedLine.getContig(),
							bedLine.getStart(),
							bedLine.getEnd()
							);
					while(r.hasNext())
						{
						final SAMRecord rec=r.next();
						if(rec.getReadUnmappedFlag()) continue;
						if(this.filter.filterOut(rec)) continue;
						if(!rec.getReferenceName().equals(bedLine.getContig())) continue;
						
						final Cigar cigar=rec.getCigar();
						if(cigar==null) continue;
			    		int refpos1=rec.getAlignmentStart();
			    		for(final CigarElement ce:cigar)
			    			{
			    			final CigarOperator op=ce.getOperator();
			    			if(!op.consumesReferenceBases()) continue;
			    			if(op.consumesReadBases())
			    				{
			    				for(int i=0;i< ce.getLength();++i)
		    		    			{
									if(refpos1+i>= bedLine.getStart() && refpos1+i<=bedLine.getEnd())
										{
										counts[refpos1+i-bedLine.getStart()]++;
										}
				    				}
			    				}
			    			refpos1+=ce.getLength();
			    			if(refpos1>bedLine.getEnd()) break;
			    			}
						}
					
					r.close();
					
					Arrays.sort(counts);
					
					int count_no_coverage=0;
					double mean=0;
					for(final int cov:counts)
						{
						if(cov<=MIN_COVERAGE) ++count_no_coverage;
						mean+=cov;
						}
					mean/=counts.length;
					
	                final double median_depth;
	                final int mid_x= counts.length/2;
	                if(counts.length%2==0)
	                        {
	                        median_depth =  (counts[mid_x-1]+counts[mid_x])/2.0;
	                        }
	                else
	                        {
	                        median_depth =  counts[mid_x];
	                        }
	
					
					pw.println(
							bedLine.getContig()+"\t"+
							(bedLine.getStart()-1)+"\t"+
							(bedLine.getEnd())+"\t"+
							counts.length+"\t"+
							(genomicSequence==null?
								"":
								String.valueOf((genomicSequence.getGCPercent(bedLine.getStart()-1,bedLine.getEnd())).getGCPercentAsInteger())+"\t")+
							counts[0]+"\t"+
							counts[counts.length-1]+"\t"+
							mean+"\t"+median_depth+"\t"+
							count_no_coverage+"\t"+
							(int)(((counts.length-count_no_coverage)/(double)counts.length)*100.0)
							);
					}
				pw.flush();
				pw.close();pw=null;
				LOG.info("done");
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
			CloserUtil.close(bedIn);
			CloserUtil.close(samReader);
			}
		}
	
	public static void main(final String[] args) throws Exception
		{
		new BamStats04().instanceMainWithExit(args);
		}

}

