package com.github.lindenb.jvarkit.util.picard;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordComparator;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloseableIterator;

/** a fast version of merging samrecorditerator, no merging of read groups
 * assuming sam sequence dictionaries are the same */
public class MergingSamRecordIterator implements Iterator<SAMRecord>,CloseableIterator<SAMRecord>
	{
	private List<PeekableIterator<SAMRecord>> iterators=new ArrayList<>();
	private SAMRecord _next=null;
	private SAMRecordComparator comparator=new SAMRecordCoordinateComparator();
	public MergingSamRecordIterator(List<SAMRecordIterator> iterators)
		{
		this.iterators=new ArrayList<PeekableIterator<SAMRecord>>(iterators.size());
		for(SAMRecordIterator it:iterators)
			{
			this.iterators.add(new PeekableIterator<SAMRecord>(it));
			}
		}
	

	@Override
	public boolean hasNext()
		{
		if(_next!=null) return true;
		if(iterators.isEmpty()) return false;
		int i=0;
		int lowest_index=-1;
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
			if(best==null || this.comparator.compare(rec, best)<0)
				{
				lowest_index=i;
				best=rec;
				}
			++i;
			}
		if(lowest_index!=-1)
			{
			_next= iterators.get(lowest_index).next();
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
		this.iterators.clear();
		this._next=null;
		}
	}
