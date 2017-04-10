package com.github.lindenb.jvarkit.util.bio.gtf;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import htsjdk.samtools.util.Locatable;

public class GTFGene implements Locatable
	{
	public class Exon implements Locatable
		{
		String exon_id;
		int start;
		int end;
		int index;
		long count_prev_and_next=0L;
		long count_prev_and_curr=0L;
		long count_curr_and_next=0L;
		long count_curr_only=0L;
		long count_others=0L;
		
		
		public GTFGene getGene()
			{
			return GTFGene.this;
			}
		@Override
		public int getStart() {
			return start;
			}
		
		@Override
		public int getEnd() {
			return end;
			}
		
		public String getContig()
			{
			return getGene().getContig();
			}
		
		public List<Exon> getPrev()
			{
			if(index==0) return Collections.emptyList();
			return getGene().exons.subList(0,index);
			}
		public List<Exon> getNext()
			{
			if(index+1>=getGene().exons.size()) return Collections.emptyList();
			return getGene().exons.subList(index+1,getGene().exons.size());
			}
		boolean contains(int pos)
			{
			return start<=pos && pos<=end;
			}
		@Override
		public String toString()
			{
			return ""+start+"-"+end;
			}
		}
	
	String chrom;
	String gene_name;
	String gene_id;
	String transcript_id;
	List<Exon> exons=new ArrayList<Exon>();
	
	GTFGene()
		{
		
		}

 	GTFGene(List<GTFLine> lines)
		{
 		for(final GTFLine item: lines) {
 			String token=item.getType();
 			if(token.equals("gene")) {
				continue;
 			}
			
			else if(token.equals("transcript")) {
				//tx = item.interval;
				continue;
				}
			
			else if(token.equals("exon")) {
				//exons.add( item.interval);
				continue;
			}
			
			else if(token.equals("CDS")) {
				//cds.add( item.interval);
				continue;
				}
 			}

		}
	
	Exon createExon(int start,int end)
		{
		Exon exon=new Exon();
		exon.start=start;
		exon.end=end;
		this.exons.add(exon);
		return exon;
		}
	
	@Override
	public String getContig() {
		return chrom;
		}
	
	@Override
	public String toString()
		{
		return transcript_id+" "+exons;
		}
	
	@Override
	public int getStart() {
		return this.exons.get(0).getStart();
		}
	
	@Override
	public int getEnd() {
		return this.exons.get(this.exons.size()-1).getEnd();
		}
	}
