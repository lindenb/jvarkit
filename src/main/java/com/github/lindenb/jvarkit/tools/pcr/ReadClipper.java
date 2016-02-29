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
package com.github.lindenb.jvarkit.tools.pcr;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.Interval;

public class ReadClipper
	{
	private static final Logger LOG=Logger.getLogger("jvarkit");
	private String programGroup = null;
	
	public void setProgramGroup(String programGroup) {
		this.programGroup = programGroup;
	}
	
	public String getProgramGroup() {
		return programGroup;
	}
	
	public SAMRecord clip(SAMRecord rec,final Interval fragment)
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

		
		int newStart = rec.getAlignmentStart();
		Cigar cigar = rec.getCigar();
		if(cigar==null)
			{
			LOG.warning("cigar missing in "+rec);
			return rec;
			}
		
		List<CigarOperator> operators = new ArrayList<>();
		//expand cigar	
		for(CigarElement ce:cigar.getCigarElements())
			{
			CigarOperator op = ce.getOperator();
			for(int x=0;x < ce.getLength();++x)
				{
				operators.add(op);
				}
			}
		
		/* 5' side */				
		if(rec.getAlignmentStart() < fragment.getStart())
			{
			int operator_index =0;
			newStart = fragment.getStart();
			int refPos1 = rec.getUnclippedStart();
			
			while(	operator_index < operators.size() &&
					refPos1< fragment.getStart())
				{
				CigarOperator op = operators.get(operator_index);
				if(op.consumesReferenceBases() )
					{
					refPos1++;
					}
				if(op.consumesReadBases())
					{
					operators.set(operator_index, CigarOperator.S);
					}
				operator_index++;
				}
			//shouln't start with a problem
			while(operator_index < operators.size())
				{
				CigarOperator op = operators.get(operator_index);
				if(op.equals(CigarOperator.M) || op.equals(CigarOperator.EQ))
					{
					break;
					}
				if(op.consumesReadBases())//insertion...
					{
					operators.set(operator_index, CigarOperator.S);
					newStart++;
					}
				operator_index++;
				}
			int x=0;
			int y=0;
			for(x=0;x< operator_index;++x)
				{
				CigarOperator op = operators.get(y);
				if(!(op.equals(CigarOperator.S) || op.equals(CigarOperator.H)))
					{
					operators.remove(y);
					}
				else
					{
					++y;
					}
				}
			
			}

		
		/* 3' side */				
		if(rec.getAlignmentEnd() > fragment.getEnd())
			{
			int operator_index = operators.size()-1;
			int refPos1 = rec.getUnclippedEnd();

			while(operator_index >=0 &&
				refPos1 > fragment.getEnd())
				{
				CigarOperator op = operators.get(operator_index);
				if(op.consumesReferenceBases() )
					{
					refPos1--;
					}
				if(op.consumesReadBases())
					{
					operators.set(operator_index, CigarOperator.S);
					}
				operator_index--;
				}
			
			//shouln't end with a problem
			while(operator_index >=0 )
				{
				CigarOperator op = operators.get(operator_index);
				if(op.equals(CigarOperator.M) || op.equals(CigarOperator.EQ))
					{
					break;
					}
				if(op.consumesReadBases())//insertion...
					{
					operators.set(operator_index, CigarOperator.S);
					}
				operator_index--;
				}
			int y=operators.size()-1;
			int len = (operators.size()-1)-operator_index;
			for(int x=0;x< len;++x)
				{
				CigarOperator op = operators.get(y);
				if(!(op.equals(CigarOperator.S) || op.equals(CigarOperator.H)))
					{
					operators.remove(y);
					}
				else
					{
					--y;
					}
				}
			
			}
		
		
		
		//build new cigar
		boolean found_M=false;
		List<CigarElement> newcigarlist=new ArrayList<>();
		int c1 = 0;
		while(c1 < operators.size())
			{
			CigarOperator op =  operators.get(c1);
			if(op.equals(CigarOperator.M) || op.equals(CigarOperator.EQ))
				{
				found_M=true;
				}
			int c2=c1;
			while(c2< operators.size() && op.equals(operators.get(c2)))
				{
				++c2;
				}
			newcigarlist.add(new CigarElement(c2-c1,op));
			c1=c2;
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
