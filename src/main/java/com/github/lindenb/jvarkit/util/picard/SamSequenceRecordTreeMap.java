/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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


*/
package com.github.lindenb.jvarkit.util.picard;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.SAMSequenceDictionary;

public class SamSequenceRecordTreeMap<T>
	{
	private SAMSequenceDictionary dict=null;
	private List<IntervalTree<T>> chroms;
	public SamSequenceRecordTreeMap(final SAMSequenceDictionary dict)
		{
		this.dict=dict;
		this.chroms=new ArrayList<IntervalTree<T>>(this.dict.size());
		}

	
	public SAMSequenceDictionary getSAMSequenceDictionary()
		{
		return dict;
		}
	
	protected IntervalTree<T> tree(int tid)
		{
		if(tid<0 || tid>=this.chroms.size()) return null;
		return this.chroms.get(tid);
		}
	
	public List<T> getOverlapping(String chrom,int start1,int end1)
		{
		return getOverlapping(getSAMSequenceDictionary().getSequenceIndex(chrom),start1,end1);
		}
	
    public List<T> getOverlapping(int tid,int start1,int end1)
    	{
        IntervalTree<T> tree=tree(tid);
        if(tree==null) return Collections.emptyList();
        List<T> L=new ArrayList<T>();
        for(Iterator<IntervalTree.Node<T>> iter = tree.overlappers(start1,end1);
        		iter.hasNext();
        		)
        		{
                L.add(iter.next().getValue());
        		}
       
        return L;
    	}
    
    public boolean containsOverlapping(String chrom,int start1,int end1)
		{
		return containsOverlapping(getSAMSequenceDictionary().getSequenceIndex(chrom),start1,end1);
		}

	public boolean  containsOverlapping(int tid,int start1,int end1)
		{
	    IntervalTree<T> tree=tree(tid);
	    return tree!=null && tree.overlappers(start1,end1).hasNext();
	 	}
	
	/** return true if there is one element defined for this chromosome */
	public boolean  containsChromosome(int tid)
		{
	    return tree(tid)!=null;
	 	}
	
	/** return true if there is one element defined for this chromosome */
	public boolean  containsChromosome(String chrom)
		{
	    return containsChromosome(getSAMSequenceDictionary().getSequenceIndex(chrom));
	 	}
	
	
	/** inserts object o at chrom/start1/end1 
	 *  returns true if object was inserted
	 * */
	public boolean put(String chrom,int start1,int end1,T o)
		{
		return put(getSAMSequenceDictionary().getSequenceIndex(chrom),start1,end1,o);
		}
	
	/** inserts object o at tid/start1/end1 
	 *  returns true if object was inserted
	 * */
	public boolean put(int tid,int start1,int end1,T o)
		{
		if(tid<0) return false;
		IntervalTree<T> m;
		if(this.chroms.size()<=tid)
			{
			if(tid>=this.getSAMSequenceDictionary().size()) return false;
			while(this.chroms.size()<tid) this.chroms.add(null);
			m=new IntervalTree<T>();
			this.chroms.add(m);
			}
		else 
			{
			m= this.chroms.get(tid);
			if(m==null)
				{
				m=new IntervalTree<T>();
				this.chroms.set(tid, m);
				}
			}
		m.put(start1, end1, o);
		return true;
		}
	
	public T get(String chrom,int start1,int end1)
		{
		return get(getSAMSequenceDictionary().getSequenceIndex(chrom),start1,end1);
		}
	
	public T get(int tid,int start1,int end1)
		{
		
		if(tid<0 || this.chroms.size()<=tid)
			{
			return null;
			}
		
		IntervalTree<T> m= this.chroms.get(tid);
		if(m==null)
			{
			return null;
			}
		IntervalTree.Node<T> node=m.find(start1, end1);
		if(node==null) return null;
		return node.getValue();
		}
	
	public boolean isEmpty()
		{
		return this.chroms.isEmpty();
		}
	}
