package com.github.lindenb.jvarkit.tools.tview;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;


import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
//import net.sf.picard.util.Log;
import net.sf.picard.util.SamLocusIterator;
import net.sf.picard.util.SamLocusIterator.RecordAndOffset;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;


public class TView
	{
	//private static final Log LOG = Log.getInstance(TView.class);
    
    private static class SamPixel
    	{
    	int refPos=-1;
    	int readPos=-1;
    	SamPixel next=null;
    	}
    
    private static class Read
    	{
    	SAMRecord record;
    	SamPixel first=null;
    	SamPixel last()
			{
			SamPixel last=first;
			while(last.next!=null)last=last.next;
			return last;
			}

    	}
    
    
	public int execute(
			final SAMFileReader samReader,
			final IndexedFastaSequenceFile reference,
			Interval interval,
			final TViewHandler handler
			)
		{   
		
		if(samReader==null) throw new NullPointerException("undefined samReader");
		if(interval==null)
			{
			SAMSequenceRecord ssr=samReader.getFileHeader().getSequence(0);
			interval=new Interval(ssr.getSequenceName(), 1, 80);
			}
		
		List<List<SAMRecord>> rows=new ArrayList<List<SAMRecord>>();
		SamLocusIterator slit=null;
		Iterator<SamLocusIterator.LocusInfo> iter=null;
		try {
	        
            //loop over reads overlaping region
            SAMRecordIterator sri=samReader.queryOverlapping(
            		interval.getSequence(),
            		interval.getStart(), 
            		interval.getEnd()
            		);
            while(sri.hasNext())
            	{
            	SAMRecord rec=sri.next();
            	if(rec.getReadUnmappedFlag()) continue;
                if(rec.getAlignmentEnd() < interval.getStart() ) continue;
                if(rec.getAlignmentStart() > interval.getEnd() ) continue;
            	
                Cigar c=rec.getCigar();
                if(c==null) continue;
                
                
                
                int readpos0=0;
                int refpos1=rec.getAlignmentStart();
                for(CigarElement ce:c.getCigarElements())
                	{
                	switch(ce.getOperator())
                		{
                		case S:
                			{
                			readpos0+=ce.getLength();
                			break;
                			}
                		case D:
                			{
                			
                			break;
                			}
                		case EQ: case M: case X:
                			{
                			readpos0+=ce.getLength();
                			refpos1+=ce.getLength();
                			break;
                			}
                		default:
                			{
                			throw new IllegalStateException("Operator not handled: "+rec.getCigarString());
                			}
                		}
                	}
            	int y=0;
                for(y=0;y< rows.size();++y)
                	{
                	List<SAMRecord> row=rows.get(y);
                	SAMRecord last=row.get(row.size()-1);
                	if(last.getAlignmentEnd()+1 < rec.getAlignmentStart())
                		{
                		row.add(rec);
                		break;
                		}
                	}
                if(y==rows.size())
                	{
                	List<SAMRecord> row=new ArrayList<SAMRecord>();
                	row.add(rec);
                	rows.add(row);
                	}
            	}
            sri.close();
            
            
            
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
	        handler.beginDocument();
	        
	        /* get REFERENCE SEQUENCE */
	        if(reference!=null)
		        {
	        	handler.beginReferenceSeq();
		        for(int refPos=interval.getStart();refPos<=interval.getEnd();++refPos)
		        	{
		        	Integer x=maxInsertAt.get(refPos);
		        	while(x!=null && x>0)
	    				{
		        		handler.reference(reference, interval.getSequence(),-1);
	    				x--;
	    				}
		        	handler.reference(reference,interval.getSequence(), refPos);
		        	}
		        handler.endReferenceSeq();
		        }
	        
	        
	        
	        for(List<Read> row:rows)
	        	{
	        	handler.beginRow();
	        	int refPos=interval.getStart();
	        	for(Read r:row)
		        	{
	        		
	        		while(refPos< r.record.getAlignmentStart())
	        			{
	        			handler.whiteSpace();
	        			Integer x=maxInsertAt.get(refPos);
	        			while(x!=null && x>0)
	        				{
	        				handler.whiteSpace();
	        				x--;
	        				}
	        			refPos++;
	        			}
	        		handler.beginSAMRecord(r.record);
		        	SamPixel pix=r.first;
		        	int last_insert_size=0;
		        	while(pix!=null)
		        		{
		        		if(pix.readPos!=-1)
		        			{
		        			if(pix.refPos==-1) last_insert_size++;
		        			}
		        		
		        		if(pix.refPos!=-1)
		        			{
		        			Integer x=maxInsertAt.get(pix.refPos);
		        			if(x!=null) x-=last_insert_size;
		        			while(x!=null && x>0 && pix!=r.first)//don't do this for first pixel
		        				{
		        				handler.deletion();
		        				x--;
		        				}
		        				
		        			refPos=pix.refPos;
		        			last_insert_size=0;
		        			}
		        		handler.base(
		        				r.record,
		        				pix.readPos,
		        				reference,
		        				pix.refPos
		        				);
		        		pix=pix.next;
		        		}
		        	handler.endSAMRecord(r.record);
		        	refPos++;//????????????????????????========
		        	}
	        	handler.endRow();
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
			}
		return 0;
		}
	}
