/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.pcr;

import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.Interval;

public class ReadClipper
	{
	private static final Logger LOG=Logger.build(ReadClipper.class).make();
	private String programGroup = null;
	public void setProgramGroup(final String programGroup) {
		this.programGroup = programGroup;
	}
	
	public String getProgramGroup() {
		return programGroup;
	}
	
	
	private static class BaseOp
		{
		CigarOperator op;
		final int refPos;
		BaseOp(final CigarOperator op,final int refPos)
			{ 
			this.op= op;
			this.refPos = refPos;
			}
		boolean isMatch()
			{
			switch(op)
				{
				case EQ: case M: return true;
				default: return false;
				}
			}		
		boolean isDeletion()
			{
			switch(op)
				{
				case D: case N: return true;
				default: return false;
				}
			}
		}
	
	public SAMRecord clip(final SAMRecord rec,final Interval fragment)
		{
		if(rec.getReadUnmappedFlag())
			{
			return rec;	
			}
		
		if(!fragment.getContig().equals(rec.getContig()))
			{
			return rec;
			}
		
		if(rec.getAlignmentEnd() < fragment.getStart())
			{
			return rec;
			}
		
		if(rec.getAlignmentStart() > fragment.getEnd())
			{
			return rec;
			}

		
		Cigar cigar = rec.getCigar();
		if(cigar==null)
			{
			LOG.warning("cigar missing in "+rec);
			return rec;
			}
		final List<BaseOp> bases = new ArrayList<>(
				cigar.getCigarElements().stream().
				mapToInt(C->C.getLength()).sum()
				);
		//expand cigar	
		int refPos= rec.getUnclippedStart();
		for(final CigarElement ce:cigar.getCigarElements())
			{
			final CigarOperator op = ce.getOperator();
			if(op.equals(CigarOperator.P)) continue;
			for(int x=0;x < ce.getLength();++x)
				{
				bases.add(new BaseOp(op,
					op.consumesReferenceBases() || op.isClipping()? refPos:-1
					));
				
				if(op.consumesReferenceBases() || op.isClipping())
					{
					refPos++;
					}
				}
			}
		
		/* 5' side */	
		int newStart = rec.getAlignmentStart();			
		int x=0;
		while(x<bases.size())
			{
			final BaseOp b = bases.get(x);
			if(b.refPos!=-1 && b.isMatch() && b.refPos>=fragment.getStart())
				{
				newStart=b.refPos;
				break;
				}
			else if(b.isDeletion())
				{
				bases.remove(x);
				continue;
				}
			else if(!b.op.isClipping())
				{
				b.op=CigarOperator.S;
				++x;
				}
			else
				{
				++x;
				}
			}
		
		/* 3' side */
		x = bases.size()-1;
		while(x>=0)
			{
			final BaseOp b = bases.get(x);
			if(b.refPos!=-1 && b.isMatch() && b.refPos<=fragment.getEnd())
				{
				break;
				}
			else if(b.isDeletion())
				{
				bases.remove(x);
				--x;
				continue;
				}
			else if(!b.op.isClipping())
				{
				b.op=CigarOperator.S;
				--x;
				}
			else
				{
				--x;
				}
			}
		
		
		//build new cigar
		boolean found_M=false;
		final List<CigarElement> newcigarlist=new ArrayList<>();
		x = 0;
		while(x < bases.size())
			{
			final CigarOperator op =  bases.get(x).op;
			if(op.equals(CigarOperator.M) || op.equals(CigarOperator.EQ))
				{
				found_M=true;
				}
			int x2=x;
			while(x2< bases.size() && op.equals(bases.get(x2).op))
				{
				++x2;
				}
			newcigarlist.add(new CigarElement(x2-x,op));
			x=x2;
			}
		
		if(!found_M)
			{
			SAMUtils.makeReadUnmappedWithOriginalTags(rec);
			if(this.programGroup!=null) {
				rec.setAttribute("PG",programGroup);
				}
			return rec;
			}
		cigar = new Cigar(newcigarlist);				
		rec.setCigar(cigar);
		rec.setAlignmentStart(newStart);
		if(this.programGroup!=null) {
			rec.setAttribute("PG",programGroup);
		}
		return rec;
		}
	}
