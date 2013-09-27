package com.github.lindenb.jvarkit.tools.bamstats04;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.regex.Pattern;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.IOUtils;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.samtools.Cigar;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

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
		SAMFileReader samReader = null;
		try
			{
			System.out.println("#chrom\tstart\tend\tlength\tmincov\tmaxcov\tmean\tnocoveragepb\tpercentcovered");
			LOG.info("Scanning "+IN);
			Pattern tab=Pattern.compile("[\t]");
			String tokens[];
			bedIn=IOUtils.openFileForBufferedReading(BEDILE);
			samReader = new SAMFileReader(IN);
			samReader.setValidationStringency(super.VALIDATION_STRINGENCY);
			String line=null;
			while((line=bedIn.readLine())!=null)
				{
				if(line.isEmpty() || line.startsWith("#")) continue;
				LOG.debug(line);
				tokens=tab.split(line,5);
				if(tokens.length<3) throw new IOException("bad bed line in "+line+" "+this.BEDILE);
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
					if(NO_ORPHAN && !rec.getProperPairFlag()) continue;
					if(!rec.getReferenceName().equals(chrom)) continue;
					if(rec.getMappingQuality()==255 && rec.getMappingQuality()< this.MMQ)
						{
						continue;
						}
					Cigar cigar=rec.getCigar();
					if(cigar==null) continue;
		    		int refpos1=rec.getUnclippedStart();
		    		for(CigarElement ce:cigar.getCigarElements())
		    			{
	    				switch(ce.getOperator())
	    					{
	    					case H:break;
	    					case S:break;
	    					case I:break;
	    					case P:break;
	    					case N:// reference skip
	    					case D://deletion in reference
	    						{
		    					refpos1+=ce.getLength();
	    						break;
	    						}
	    					case M:
	    					case EQ:
	    					case X:
	    						{
	    						for(int i=0;i< ce.getLength() && refpos1<= chromEnd1;++i)
		    		    			{
	    							if(refpos1>= chromStart1 && refpos1<=chromEnd1)
	    								{
	    								counts[refpos1-chromStart1]++;
	    								}
		    						refpos1++;
	    		    				}
	    						break;
	    						}
	    					default: throw new IllegalStateException(
	    							"Doesn't know how to handle cigar operator:"+ce.getOperator()+
	    							" cigar:"+cigar
	    							);

	    					}
		    				
		    			}
					}
				
				r.close();
				
				Arrays.sort(counts);
				
				int count_no_coverage=0;
				double mean=0;
				for(int cov:counts)
					{
					if(cov<=MIN_COVERAGE) ++count_no_coverage;
					mean+=cov;
					}
				mean/=counts.length;
				
				System.out.println(
						tokens[0]+"\t"+tokens[1]+"\t"+tokens[2]+"\t"+
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
		if(bedIn!=null) try {bedIn.close();}catch(IOException err) {}
		if(samReader!=null) samReader.close();
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

