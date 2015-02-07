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
package com.github.lindenb.jvarkit.util.picard;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.util.CloseableIterator;

/** a fast ? version of merging samrecorditerator, no merging of read groups
 * assuming sam sequence dictionaries are the same */
public class MergingSamRecordIterator implements Iterator<SAMRecord>,CloseableIterator<SAMRecord>
	{
	private List<PeekableIterator<SAMRecord>> iterators=new ArrayList<PeekableIterator<SAMRecord>>();
	private SAMRecord _next=null;
	private Comparator<SAMRecord> comparator;
	private ArrayList<Integer> bestIndexes;
	public MergingSamRecordIterator(
			Comparator<SAMRecord> comparator,
			List<CloseableIterator<SAMRecord>> iterators
			)
		{
		this.iterators=new ArrayList<PeekableIterator<SAMRecord>>(iterators.size());
		for(CloseableIterator<SAMRecord> it:iterators)
			{
			this.iterators.add(new PeekableIterator<SAMRecord>(it));
			}
		this.bestIndexes=new ArrayList<>();
		this.comparator=comparator;
		if(this.comparator==null) throw new NullPointerException("comparator is null");
		}
	public MergingSamRecordIterator(List<CloseableIterator<SAMRecord>> iterators)
		{
		this(new SAMRecordCoordinateComparator(),iterators);
		}
	

	@Override
	public boolean hasNext()
		{
		if(_next!=null) return true;
		if(!this.bestIndexes.isEmpty())
			{
			int las_index= this.bestIndexes.remove(this.bestIndexes.size()-1);
			_next= iterators.get(las_index).next();
			return true;
			}	
		if(iterators.isEmpty()) return false;
		int i=0;
		SAMRecord best=null;
		while(i<iterators.size())
			{
			PeekableIterator<SAMRecord> peeker=iterators.get(i);
			if(!peeker.hasNext())
				{
				peeker.close();
				iterators.remove(i);
				continue;
				}
			SAMRecord rec=peeker.peek();
		
			
			if(best==null)
				{
				this.bestIndexes.clear();
				this.bestIndexes.add(i);
				best=rec;
				}
			else
				{
				int cmp=this.comparator.compare(rec, best);
				if(cmp<0)
					{
					this.bestIndexes.clear();
					this.bestIndexes.add(i);
					best=rec;
					}
				else if(cmp==0)
					{
					this.bestIndexes.add(0,i);//push front
					}
				}
			++i;
			}
		if(!bestIndexes.isEmpty())
			{
			int las_index= this.bestIndexes.remove(this.bestIndexes.size()-1);
			_next= iterators.get(las_index).next();
			}
		return _next!=null;
		}

	@Override
	public SAMRecord next()
		{
		if(!hasNext()) throw new NoSuchElementException();
		SAMRecord o=_next;
		_next=null;
		return o;
		}

	@Override
	public void remove()
		{
		throw new UnsupportedOperationException();
		}
	@Override
	public void close()
		{
		for(PeekableIterator<SAMRecord> it:this.iterators)
			{
			it.close();
			}
		this.bestIndexes.clear();
		this.iterators.clear();
		this._next=null;
		}
	}
