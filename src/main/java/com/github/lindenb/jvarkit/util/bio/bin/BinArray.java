package com.github.lindenb.jvarkit.util.bio.bin;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


public class BinArray<T> implements Iterable<T>
	{
	private static class Node<T>
		{
		T value=null;
		int start0;
		int end0;
		Node<T> next=null;
		Node(int start0,int end0,T value)
			{
			this.start0=start0;
			this.end0=end0;
			this.value=value;
			}
		}
	
	private List<Node<T>> bin2node=new ArrayList<Node<T>>();

	public boolean put(int start0,int end0,T value)
		{
		Node<T> node=new Node<T>(start0,end0,value);
		int index=binFromRangeStandard(start0,end0);
		while(index>=bin2node.size()) bin2node.add(null);
		Node<T> root=bin2node.get(index);
		if(root==null)
			{
			bin2node.set(index,node);
			}
		else
			{
			while(root.next!=null)
				{
				root=root.next;
				}
			root.next=node;
			}
		return true;
		}
	
	 public boolean containsOverlapping(int start0,int end0)
	 	{
		return overlapping(start0,end0).hasNext();
	 	}
	
	 public Iterator<T> overlapping(int start0,int end0)
	 	{
		return new Iter(start0, end0);
	 	}
	 
	 
	 public List<T> getOverlapping(int start0,int end0)
	 	{
		List<T> L=new ArrayList<T>();
		Iterator<T> r=overlapping(start0, end0);
		while(r.hasNext())
			{
			L.add(r.next());
			}
		return L;
	 	}
	 
	 @Override
	public Iterator<T> iterator() {
		return overlapping(0,Integer.MAX_VALUE);
	 	}
	 
	 private class Iter
	 	implements Iterator<T>
	 	{
		private boolean _has_next_called=false;
		private boolean _has_next=false;
		private T _next=null;
		private int beg;
		private int end;
		private BinIterator binIterator;
		private Node<T> currNode=null;
		private int binIndex=-1;
		
		Iter(int beg,int end)
			{
			this.beg=beg;
			this.end=end;
			this.binIterator=new BinIterator(beg, end);
			}		
		protected boolean overlap(Node<T> n)
			{
			if(n.start0>=this.end) return false;
			if(n.end0<=this.beg) return false;
			return true;
			}
		private boolean shift()
			{
			for(;;)
				{
				if(currNode!=null)
					{
					Node<T> tmp=currNode;
					currNode=currNode.next;
					
					if(overlap(tmp))
						{
						_next=tmp.value;
						return true;
						}
					
					}
				else if(binIndex==-1)
					{
					if(!binIterator.hasNext())
						{
						return false;
						}	
					binIndex=binIterator.next();
					}
				else if(binIndex>=bin2node.size())
					{
					return false;
					}
				else
					{
					if(binIterator.hasNext())
						{
						binIndex=binIterator.next();
						if(binIndex< bin2node.size())
							{
							currNode=bin2node.get(binIndex);
							}
						}	
					
					}
				}
			}
		
	 	public boolean hasNext()
		 	{
			if(_has_next_called) return _has_next;
			_has_next_called=true;
			_has_next=shift();
			return _has_next;
		 	}
		public T next()
			{
			if(!_has_next_called) hasNext();
			if(!_has_next) throw new IllegalStateException();
			T tmp=_next;
			_next=null;
			_has_next_called=false;
			_has_next=false;
			return tmp;
			}
		@Override
		public void remove() {
			throw new UnsupportedOperationException();
			}
	 	}
	
	
	
	

static final int binOffsets[] = {512+64+8+1, 64+8+1, 8+1, 1, 0};
static final int  _binFirstShift=17;       /* How much to shift to get to finest bin. */
static final int _binNextShift=3 ;       /* How much to shift to get to next larger bin.


/* Copied from http://genomewiki.ucsc.edu/index.php/Bin_indexing_system
 * Given start,end in chromosome coordinates assign it
 * a bin.   There's a bin for each 128k segment, for each
 * 1M segment, for each 8M segment, for each 64M segment,
 * and for each chromosome (which is assumed to be less than
 * 512M.)  A range goes into the smallest bin it will fit in. */
static private int binFromRangeStandard(int start, int end)
	{
	int startBin = start, endBin = end-1;
	startBin >>= _binFirstShift;
	endBin >>= _binFirstShift;
	for (int binOffset : binOffsets)
	    {
	    if (startBin == endBin)
	        return binOffset + startBin;
	    startBin >>= _binNextShift;
	    endBin >>= _binNextShift;
	    }
	throw new IllegalStateException("start "+start+", end "+end+" out of range in findBin (max is 512M)");
	}

}
