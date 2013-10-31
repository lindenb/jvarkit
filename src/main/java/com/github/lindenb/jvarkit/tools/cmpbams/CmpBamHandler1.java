package com.github.lindenb.jvarkit.tools.cmpbams;

import java.io.PrintStream;

import com.github.lindenb.jvarkit.tools.cmpbams.entities.BamRecord;
import com.github.lindenb.jvarkit.tools.cmpbams.entities.Records;
import com.github.lindenb.jvarkit.tools.cmpbams.entities.Comparebams;
import com.github.lindenb.jvarkit.util.Counter;

public class CmpBamHandler1 extends AbstractCompareBamsHandler
	{
	private Counter<String> counter=new Counter<String>();
	@Override
	protected boolean handleRecords(PrintStream out, Records records)
		{
		
		BamRecord r1=null;
		BamRecord r2=null;
		for(BamRecord r:records.getRecord())
			{
			if(r.getFileIndex()==1)
				{
				r1=r;
				}
			else if(r.getFileIndex()==2)
				{
				r2=r;
				}
			}
		
		this.counter.incr("total");
		
		if(r1==null && r2!=null)
			{
			this.counter.incr("in BAM2, not in BAM1");
			System.err.println("MISSING\t"+records.getName());
			return true;
			}	
		else if(r2==null && r1!=null)
			{
			this.counter.incr("in BAM1, not in BAM2");
			System.err.println("MISSING\t"+records.getName());
			return true;
			}
		
		if(!r1.getFlag().equals(r2.getFlag()))
			{
			this.counter.incr("Flag change from "+r1.getFlag()+" to "+r2.getFlag());
			}
		
		if(r1.getChrom()==null && r2.getChrom()==null)
			{
			this.counter.incr("Both unmapped");
			return true;

			}
		
		if( (r1.getChrom()!=null && !r1.getChrom().equals(r2.getChrom())) ||
			(r1.getChrom()==null && r2.getChrom()!=null)
			)
			{
			this.counter.incr("From BAM1:chrom "+r1.getChrom() +" to BAM2:"+r2.getChrom() );
			return true;
			}
		
		//from here both mapped
		if(!r1.getPos().equals(r2.getPos()))
			{
			int length=Math.abs(r1.getPos()-r2.getPos());
			this.counter.incr("Same Chrom but Position-shift: "+ ((int)(length/100.0))*10);
			}
		
		else if(!r1.getCigar().equals(r2.getCigar()))
			{
			this.counter.incr("Same Position. Cigar changed");
			}
		else
			{
			this.counter.incr("No Change");
			}
		
		return true;
		}
	
	@Override
	protected void setHeader(Comparebams.Header header) {
		super.setHeader(header);		
		}
	
	@Override
	protected void finish(PrintStream out) {
		for(String s:counter.keySet())
			{
				out.print(s);
				out.print('\t');
				out.print(counter.count(s));
				out.print('\t');
				out.print(String.format("%.2f%%",((counter.count(s)/(double)counter.count("total"))*100.0)));
				out.println();
			}
		}
	
	
		/**
	 * @param args
	 */
	public static void main(String[] args) {
		new CmpBamHandler1().instanceMainWithExit(args);

	}

}
