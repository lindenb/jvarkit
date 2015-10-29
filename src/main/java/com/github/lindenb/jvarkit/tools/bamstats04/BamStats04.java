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
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collection;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class BamStats04 extends AbstractBamStats04
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(BamStats04.class);

	@Override
	public Command createCommand()
		{
		return new MyCommand();
		}
	
	static private class MyCommand extends AbstractBamStats04.AbstractBamStats04Command
		{    
		
		@SuppressWarnings("resource")
		@Override
		protected Collection<Throwable> call(String IN) throws Exception
			{			
			if( BEDILE==null) return wrapException("undefined bed file");
			PrintStream pw = null;
			BufferedReader bedIn=null;
			SamReader samReader = null;
			try
				{
				pw =  openFileOrStdoutAsPrintStream();
				pw.println("#chrom\tstart\tend\tlength\tmincov\tmaxcov\tmean\tnocoveragepb\tpercentcovered");
				LOG.info("Scanning "+IN);
				Pattern tab=Pattern.compile("[\t]");
				String tokens[];
				bedIn=IOUtils.openFileForBufferedReading(BEDILE);
				
				samReader = openSamReader(IN);
				String line=null;
				while((line=bedIn.readLine())!=null)
					{
					if(line.isEmpty() || line.startsWith("#")) continue;
					LOG.debug(line);
					tokens=tab.split(line,5);
					if(tokens.length<3)
						{
						samReader.close();
						bedIn.close();
						return wrapException("bad bed line in "+line+" "+this.BEDILE);
						}
					String chrom=tokens[0];
					int chromStart1= 1+Integer.parseInt(tokens[1]);//add one
					int chromEnd1= Integer.parseInt(tokens[2]);
					/* picard javadoc:  - Sequence name - Start position (1-based) - End position (1-based, end inclusive)  */
					
					
					int counts[]=new int[chromEnd1-chromStart1+1];
					if(counts.length==0) continue;
					Arrays.fill(counts, 0);
					
					/**
					 *     start - 1-based, inclusive start of interval of interest. Zero implies start of the reference sequence.
		    		*	   end - 1-based, inclusive end of interval of interest. Zero implies end of the reference sequence. 
					 */
					SAMRecordIterator r=samReader.queryOverlapping(chrom, chromStart1, chromEnd1);
					while(r.hasNext())
						{
						SAMRecord rec=r.next();
						if(rec.getReadUnmappedFlag()) continue;
						if(NO_VENDOR && rec.getReadFailsVendorQualityCheckFlag()) continue;
						if(NO_DUP && rec.getDuplicateReadFlag() ) continue;
						if(rec.getReadPairedFlag())
							{
							if(NO_ORPHAN && !rec.getProperPairFlag()) continue;
							}
						if(rec.isSecondaryOrSupplementary()) continue;
						if(!rec.getReferenceName().equals(chrom)) continue;
						if(rec.getMappingQuality()==255 ||
							rec.getMappingQuality()==0 ||
							rec.getMappingQuality()< super.MIN_MAPING_QUALITY)
							{
							continue;
							}
						Cigar cigar=rec.getCigar();
						if(cigar==null) continue;
			    		int refpos1=rec.getAlignmentStart();
			    		for(CigarElement ce:cigar.getCigarElements())
			    			{
			    			CigarOperator op=ce.getOperator();
			    			if(!op.consumesReferenceBases()) continue;
			    			if(op.consumesReadBases())
			    				{
			    				for(int i=0;i< ce.getLength();++i)
		    		    			{
									if(refpos1+i>= chromStart1 && refpos1+i<=chromEnd1)
										{
										counts[refpos1+i-chromStart1]++;
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
					for(int cov:counts)
						{
						if(cov<=super.MIN_COV) ++count_no_coverage;
						mean+=cov;
						}
					mean/=counts.length;
					
					pw.println(
							tokens[0]+"\t"+tokens[1]+"\t"+tokens[2]+"\t"+
							counts.length+"\t"+
							counts[0]+"\t"+
							counts[counts.length-1]+"\t"+
							mean+"\t"+
							count_no_coverage+"\t"+
							(int)(((counts.length-count_no_coverage)/(double)counts.length)*100.0)
							);
					}
				return RETURN_OK;
				}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(pw);
			CloserUtil.close(bedIn);
			CloserUtil.close(samReader);
			}
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

