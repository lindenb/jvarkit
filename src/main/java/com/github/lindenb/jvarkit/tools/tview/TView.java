package com.github.lindenb.jvarkit.tools.tview;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.github.lindenb.jvarkit.util.picard.IntervalUtils;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.picard.util.Log;
import net.sf.picard.util.SamLocusIterator;
import net.sf.picard.util.SamLocusIterator.RecordAndOffset;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;


public class TView extends CommandLineProgram
	{
	private static final Log LOG = Log.getInstance(TView.class);
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Sequence logo for different alleles or generated from SAM/BAM http://www.biostars.org/p/73021";
    @Option(shortName= StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Genome Reference",optional=false)
    public File REF=null;

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="A BAM file to process.")
    public File INPUT=null;
    @Option(shortName= "L", doc="Region to observe: chrom:start-end",optional=false)
    public String REGION="";
    
    private static class SamPixel
    	{
    	int refPos=-1;
    	int readPos=-1;
    	SamPixel next=null;
    	}
    
    private static class Read
    	{
    	int x=0;
    	SAMRecord record;
    	SamPixel first=null;
    	SamPixel last()
			{
			SamPixel last=first;
			while(last.next!=null)last=last.next;
			return last;
			}
    	

    	}
    
    private class Handler
    	{
    	public void beginRow() {	 }
    	public void endRow() {	 }
    	public void beginDocument() {	 }
    	public void endDocument() {	 }
    	public void beginSAMRecord(final SAMRecord rec) {	 }
    	public void endSAMRecord(final SAMRecord rec) {	 }
    	public void whiteSpace() {}
    	public void deletion() {}
    	public void base(SAMRecord rec,int refPos,int readPos){}
    	}
    private class TTHandler extends Handler
    	{
    	public void beginRow() {	 }
    	public void endRow() {	System.out.println(); }
    	public void whiteSpace() { System.out.print(' ');}
    	public void deletion() { System.out.print('*');}
    	@Override
    	public void endDocument() {
    		System.out.println();
    		}
    	public void base(SAMRecord rec,int refPos,int readPos)
    		{	
    		System.out.print((char)rec.getReadBases()[readPos]);
    		}
    	}
    
    @Override
	protected int doWork()
		{
    	
		PrintWriter out=new PrintWriter(System.out);
		SAMFileReader samReader=null;
		SamLocusIterator slit=null;
		Iterator<SamLocusIterator.LocusInfo> iter=null;
		try {
	        samReader=new SAMFileReader(INPUT);
	        samReader.setValidationStringency(super.VALIDATION_STRINGENCY);
	        
			Interval interval=IntervalUtils.parseOne(samReader.getFileHeader().getSequenceDictionary(),REGION);
			if(interval==null)
				{
				LOG.error("Bad interval "+interval);
				return -1;
				}

	        
	        Map<SAMRecord, Read> record2seq =new LinkedHashMap<SAMRecord, Read>();
            IntervalList  iL=new  IntervalList(samReader.getFileHeader());
            iL.add(interval);
	        slit=new  SamLocusIterator(samReader,iL,true);
	        slit.setEmitUncoveredLoci(false);
	        Map<Integer,Integer> maxInsertAt=new HashMap<Integer, Integer>();

	        for(iter=slit.iterator();
                    iter.hasNext();
                    )
                {
                SamLocusIterator.LocusInfo  locusInfo=iter.next();
                int biggest_insert_here=0;
                
                if(locusInfo.getPosition() < interval.getStart() ) continue;
                if(locusInfo.getPosition() > interval.getEnd() ) continue;
                for(RecordAndOffset rao: locusInfo.getRecordAndPositions())
                	{
                	SAMRecord rec=rao.getRecord();

                	Read b=record2seq.get(rec);
                	
                	if(b==null)
                		{
                		b=new Read();
                		b.record=rec;
                		record2seq.put(rec, b);
                		}
                	SamPixel pixel=new SamPixel();
                	pixel.refPos=locusInfo.getPosition();
                	pixel.readPos=rao.getOffset();
                	if(b.first==null)
                		{
                		b.first=pixel;
                		}
                	else
                		{
                		SamPixel last=b.last();
                		//create inserts in SamRecord
                		int insert_size=0;
                		while(last.readPos +1 < pixel.readPos)
                			{
                			SamPixel insert=new SamPixel();
                			insert.readPos=last.readPos+1;
                			last.next=insert;
                			last=insert;
                			++insert_size;
                			}
                		biggest_insert_here=Math.max(insert_size, biggest_insert_here);
                		last.next=pixel;
                		}
                	
                	}
	           if(biggest_insert_here>0)
	               	{
	        	    maxInsertAt.put(locusInfo.getPosition(),biggest_insert_here);
	               	}
                
                }
	        slit.close();
	        slit=null;
	        
	        //pack
	        List<List<Read>> rows=new ArrayList<List<Read>>();
	        for(Read r:record2seq.values())
	        	{
	        	int i=0;
	        	for(i=0;i< rows.size();++i)
	        		{
	        		List<Read> row=rows.get(i);
	        		Read last=row.get(row.size()-1);
	        		if(last.record.getAlignmentEnd()+1 < r.record.getAlignmentStart())
	        			{
	        			row.add(r);
	        			break;
	        			}
	        		}
	        	if(i==rows.size())
	        		{
	        		List<Read> row=new ArrayList<TView.Read>();
	        		row.add(r);
	        		rows.add(row);
	        		}
	        	}
	        Handler handler=new TTHandler();
	        handler.beginDocument();
	        
	        for(List<Read> row:rows)
	        	{
	        	handler.beginRow();
	        	int refPos=interval.getStart();
	        	for(Read r:row)
		        	{
	        		
	        		while(refPos< r.record.getAlignmentStart())
	        			{
	        			handler.whiteSpace();
	        			//System.out.print(" ");
	        			Integer x=maxInsertAt.get(refPos);
	        			while(x!=null && x>0)
	        				{
	        				handler.whiteSpace();
	        				//System.out.print(".");
	        				x--;
	        				}
	        			refPos++;
	        			}
	        		handler.beginSAMRecord(r.record);
		        	SamPixel pix=r.first;
		        	int last_insert_size=0;
		        	while(pix!=null)
		        		{
		        		char c='*';
		        		if(pix.readPos!=-1)
		        			{
		        			c=(char)r.record.getReadBases()[pix.readPos];
		        			if(pix.refPos==-1) last_insert_size++;
		        			}
		        		
		        		if(pix.refPos!=-1)
		        			{
		        			Integer x=maxInsertAt.get(pix.refPos);
		        			if(x!=null) x-=last_insert_size;
		        			while(x!=null && x>0 && pix!=r.first)//don't do this for first pixel
		        				{
		        				handler.deletion();
		        				//System.out.print("-");
		        				x--;
		        				}
		        				
		        			refPos=pix.refPos;
		        			last_insert_size=0;
		        			}
		        		//System.out.print(c);
		        		handler.base(r.record, pix.refPos, pix.readPos);
		        		pix=pix.next;
		        		}
		        	handler.endSAMRecord(r.record);
		        	}
	        	handler.endRow();
	        	//System.out.println();
	        	}
	        handler.endDocument();
			} 
		catch (Exception e) {
			e.printStackTrace();
			return -1;
			}
		finally
			{
			if(slit!=null) slit.close();
			if(samReader!=null) samReader.close();
			out.flush();
			}
		return 0;
		}

public static void main(String[] args)
	{
	args=new String[]{"I=/home/lindenb/package/samtools-0.1.18/examples/sorted.bam","L=seq2:750-800"};
	new TView().instanceMainWithExit(args);
	}
}
