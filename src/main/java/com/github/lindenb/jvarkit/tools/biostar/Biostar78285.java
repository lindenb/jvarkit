package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.util.Arrays;
import java.util.Iterator;


import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.picard.util.SamLocusIterator;
import net.sf.picard.util.SamLocusIterator.LocusInfo;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;

public class Biostar78285 extends AbstractCommandLineProgram
	{
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+
		" Extract regions of genome that have 0 coverage See http://www.biostars.org/p/78285/ .";
	private static final Log LOG=Log.getInstance(Biostar78285.class);



    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME,doc="BAM file (sorted on coordinate). Default:stdin",optional=true)
    public File IN=null;
    @Option(shortName="SLI",doc="use a Sam locus  net.sf.picard.util.SamLocusIterator: slower but will scan the CIGAR string & detect the gaps in the reads.",optional=true)
    public boolean USE_SAMLOCUSITERATOR=false;
    
    @Override
    protected int doWork()
    	{
    	SAMFileReader samFileReader=null;
    	SAMRecordIterator iter=null;
    	SamLocusIterator sli=null;
    	try
	    	{
	    	
	    	
	    	if(IN==null)
	    		{
	    		LOG.info("read stdin");
	    		samFileReader=new SAMFileReader(System.in,false);
	    		}
	    	else
	    		{
	    		LOG.info("read "+IN);
	    		samFileReader=new SAMFileReader(IN);
	    		}
	    	samFileReader.setValidationStringency(super.VALIDATION_STRINGENCY);
	    	SAMFileHeader header=samFileReader.getFileHeader();
	    	if(header.getSortOrder()!=SortOrder.coordinate)
	    		{
	    		switch(super.VALIDATION_STRINGENCY)
	    			{
	    			case STRICT:LOG.error("SORT ORDER IS NOT 'coordinate':"+header.getSortOrder()+" (VALIDATION_STRINGENCY is STRICT)");return -1;
	    			case SILENT:break;
	    			default:LOG.warn("SORT ORDER IS NOT 'coordinate':"+header.getSortOrder());break;
	    			}
	    		}
	    	SAMSequenceDictionary dict=header.getSequenceDictionary();
	    	if(dict==null)
	    		{
	    		System.err.println("SamFile doesn't contain a SAMSequenceDictionary.");
	    		return -1;
	    		}
	    	
	    	boolean seen_tid[]=new boolean[dict.getSequences().size()];
	    	Arrays.fill(seen_tid, false);
	    	
	    	
	    	
	    	
	    	if(USE_SAMLOCUSITERATOR)
	    		{
	    		sli=new SamLocusIterator(samFileReader);
	    		sli.setEmitUncoveredLoci(true);
	    		Iterator<LocusInfo> liter=sli.iterator();
	    		LocusInfo locusInfo1=null;
	    		for(;;)
	    			{
	    			if(locusInfo1==null)
	    				{
	    				if(!liter.hasNext()) break;
	    				locusInfo1=liter.next();
	    				}
	    			int tid=locusInfo1.getSequenceIndex();
	    			//System.out.println(locusInfo1.getSequenceName()+" "+locusInfo1.getPosition()+" "+locusInfo1.getRecordAndPositions().size());
	    			seen_tid[tid]=true;
	    			if(!locusInfo1.getRecordAndPositions().isEmpty()) 
	    				{
	    				locusInfo1=null;
	    				continue;
	    				}
	    			
	    			String seqName1=locusInfo1.getSequenceName();
    				int gap_start1=locusInfo1.getPosition();
    				LocusInfo locusInfo2=null;
    				if(!liter.hasNext())
    					{
						SAMSequenceRecord ssr=dict.getSequence(tid);
		    			System.out.println(seqName1+"\t"+(gap_start1-1)+"\t"+ssr.getSequenceLength());
    					}
    				else
	    				{
	    				while(liter.hasNext())
	    					{
	    					locusInfo2=liter.next();
	    					if(locusInfo2.getSequenceIndex()!=tid)
	    						{
	    						SAMSequenceRecord ssr=dict.getSequence(tid);
	    		    			System.out.println(seqName1+"\t"+(gap_start1-1)+"\t"+ssr.getSequenceLength());
	    		    			break;
	    						}
	    	    			if(!locusInfo2.getRecordAndPositions().isEmpty())
	    	    				{
	    		    			System.out.println(seqName1+"\t"+(gap_start1-1)+"\t"+(locusInfo2.getPosition()-1));
	    		    			locusInfo2=null;
	    		    			break;
	    	    				}
	
		    				}
	    				}
    				locusInfo1=locusInfo2;
	    			
	    			}
	    	
	    		}
	    	else
		    	{
	    		int prev_tid=-1;
		    	int prev_pos1=1;
		    	
		    	iter=samFileReader.iterator();
		    	while(iter.hasNext())
		    		{
		    		SAMRecord rec=iter.next();
		    		if(rec.getReadUnmappedFlag()) continue;
		    		int tid=rec.getReferenceIndex();
		    		if(prev_tid!=-1 && prev_tid!=tid) /* chromosome has changed */
		    			{
		    			SAMSequenceRecord ssr=dict.getSequence(prev_tid);
		    			if(prev_pos1-1 < ssr.getSequenceLength())
							{
			    			System.out.println(ssr.getSequenceName()+"\t"+(prev_pos1-1)+"\t"+ssr.getSequenceLength());
			    			}
		    			prev_pos1=1;
		    			}
		    		if(prev_pos1< rec.getAlignmentStart()) /* there is a gap */
						{
		    			System.out.println(rec.getReferenceName()+"\t"+(prev_pos1-1)+"\t"+rec.getAlignmentStart());
						}
		    		seen_tid[tid]=true;
		    		prev_tid=tid;
		    		prev_pos1=Math.max(prev_pos1,rec.getAlignmentEnd()+1);
		    		}
		    	
		    	/* last reference */
				if(prev_tid!=-1 )
					{
					SAMSequenceRecord ssr=dict.getSequence(prev_tid);
					if(prev_pos1-1 < ssr.getSequenceLength())
						{
		    			System.out.println(ssr.getSequenceName()+"\t"+(prev_pos1-1)+"\t"+ssr.getSequenceLength());
						}
					}
		    	
		    	}
	    	
				
				
	    	/* unseen chromosomes */
	    	for(int i=0;i< seen_tid.length;++i)
	    		{
	    		if(seen_tid[i]) continue;
	    		SAMSequenceRecord ssr=dict.getSequence(i);
    			System.out.println(ssr.getSequenceName()+"\t0\t"+ssr.getSequenceLength());
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
    		if(iter!=null) iter.close();
    		if(sli!=null) sli.close();
    		if(samFileReader!=null) samFileReader.close();
    		}
    	
    	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Biostar78285().instanceMainWithExit(args);

	}

}
