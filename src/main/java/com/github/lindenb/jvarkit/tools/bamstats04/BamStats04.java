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
import java.io.IOException;
import java.util.Arrays;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.cmdline.Option;
import com.github.lindenb.jvarkit.util.picard.cmdline.StandardOptionDefinitions;
import com.github.lindenb.jvarkit.util.picard.cmdline.Usage;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamStats04 extends AbstractCommandLineProgram
	{
	private static final Log LOG=Log.getInstance(BamStats04.class);
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" Coverage statistics for a BED file. It uses the Cigar string instead of the start/end to get the voverage";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="BAM file to process.",
    		optional=false)
	public File IN=null;

    @Option(shortName= "BED", doc="BED File.",
    		optional=false)
	public File BEDILE=null;

    @Option(shortName= "NODUP", doc="discard duplicates", optional=false)
    public boolean NO_DUP=true;
    
    @Option(shortName= "NOORPHAN", doc="discard not properly paired", optional=false)
	public boolean NO_ORPHAN=true;
    
    @Option(shortName= "NOVENDOR", doc="discard failing Vendor Quality", optional=false)
	public boolean NO_VENDOR=true;
    
    @Option(shortName= "MIN_MAPING_QUALITY", doc="min mapping quality", optional=false)
	public int MMQ=0;

    @Option(shortName= "MIN_COV", doc="min coverage to say the position is not covered", optional=false)
	public int MIN_COVERAGE=0;

    
    ///private boolean skipDuplicates=true;
	//private int minQual=0;
	//private int basesperbin=10;
	//private int num_bin=20;
	//private boolean cumulative=true;
	
	@Override
	protected int doWork()
		{
		BufferedReader bedIn=null;
		SamReader samReader = null;
		try
			{
			System.out.println("#chrom\tstart\tend\tlength\tmincov\tmaxcov\tmean\tnocoveragepb\tpercentcovered");
			LOG.info("Scanning "+IN);
			bedIn=IOUtils.openFileForBufferedReading(BEDILE);
			
			final BedLineCodec codec= new BedLineCodec();
			SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(super.VALIDATION_STRINGENCY);
			samReader = srf.open(IN);
			String line=null;
			while((line=bedIn.readLine())!=null)
				{
				if(line.isEmpty() || line.startsWith("#")) continue;
				final BedLine bedLine = codec.decode(line);
				if(bedLine==null) continue;
				
				/* picard javadoc:  - Sequence name - Start position (1-based) - End position (1-based, end inclusive)  */
				
				
				int counts[]=new int[bedLine.getEnd()-bedLine.getStart()+1];
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
					if(NO_VENDOR && rec.getReadFailsVendorQualityCheckFlag()) continue;
					if(NO_DUP && rec.getDuplicateReadFlag() ) continue;
					if(rec.getReadPairedFlag())
						{
						if(NO_ORPHAN && !rec.getProperPairFlag()) continue;
						}
					if(rec.isSecondaryOrSupplementary()) continue;
					if(!rec.getReferenceName().equals(bedLine.getContig())) continue;
					if(rec.getMappingQuality()==255 ||
						rec.getMappingQuality()==0 ||
						rec.getMappingQuality()< this.MMQ)
						{
						continue;
						}
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
				
				System.out.println(
						bedLine.getContig()+"\t"+
						(bedLine.getStart()-1)+"\t"+
						(bedLine.getEnd())+"\t"+
						counts.length+"\t"+
						counts[0]+"\t"+
						counts[counts.length-1]+"\t"+
						mean+"\t"+
						count_no_coverage+"\t"+
						(int)(((counts.length-count_no_coverage)/(double)counts.length)*100.0)
						);
				}
			return 0;
			}
	catch(Exception err)
		{
		LOG.error(err);
		return -1;
		}
	finally
		{
		CloserUtil.close(bedIn);
		CloserUtil.close(samReader);
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception
		{
		new BamStats04().instanceMain(args);
		}

}

