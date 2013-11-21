package com.github.lindenb.jvarkit.util.bio.bin;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import net.sf.samtools.SAMSequenceDictionary;

public class SamSequenceRecordBinMap<T>
	{
	private SAMSequenceDictionary dict=null;
	private List<BinArray<T>> chroms;
	public SamSequenceRecordBinMap(SAMSequenceDictionary dict)
		{
		this.dict=dict;
		this.chroms=new ArrayList<BinArray<T>>(this.dict.size());
		}

	
	public SAMSequenceDictionary getSAMSequenceDictionary()
		{
		return dict;
		}
	
	protected BinArray<T> tree(int tid)
		{
		if(tid<0 || tid>=this.chroms.size()) return null;
		return this.chroms.get(tid);
		}
	
	public Iterator<T> overlapping(String chrom,int start0,int end0)
		{
		return overlapping(
				getSAMSequenceDictionary().getSequenceIndex(chrom),
				start0,
				end0
				);
		}
	
	public Iterator<T> overlapping(int tid,int start0,int end0)
		{
		BinArray<T> tree=tree(tid);
		if(tree==null) return new ArrayList<T>().iterator();
	    return tree.overlapping(start0, end0);
		}

	
	public List<T> getOverlapping(String chrom,int start0,int end0)
		{
		return getOverlapping(
				getSAMSequenceDictionary().getSequenceIndex(chrom),
				start0,
				end0
				);
		}
	
    public List<T> getOverlapping(int tid,int start0,int end0)
    	{
    	List<T> L=new ArrayList<T>();
        
        for(Iterator<T> iter = overlapping(tid,start0, end0);
        		iter.hasNext();
        		)
        		{
                L.add(iter.next());
        		}
       
        return L;
    	}
    
    public boolean containsOverlapping(String chrom,int start0,int end0)
		{
		return containsOverlapping(getSAMSequenceDictionary().getSequenceIndex(chrom),start0,end0);
		}

	public boolean  containsOverlapping(int tid,int start0,int end0)
		{
	    return overlapping(tid, start0, end0).hasNext();
	 	}
	
	/** inserts object o at chrom/start1/end1 
	 *  returns true if object was inserted
	 * */
	public boolean put(String chrom,int start0,int end0,T o)
		{
		return put(getSAMSequenceDictionary().getSequenceIndex(chrom),start0,end0,o);
		}
	
	/** inserts object o at tid/start1/end1 
	 *  returns true if object was inserted
	 * */
	public boolean put(int tid,int start0,int end0,T o)
		{
		if(tid<0) return false;
		BinArray<T> m;
		if(this.chroms.size()<=tid)
			{
			if(tid>=this.getSAMSequenceDictionary().size()) return false;
			while(this.chroms.size()<tid) this.chroms.add(null);
			m=new BinArray<T>();
			this.chroms.add(m);
			}
		else 
			{
			m= this.chroms.get(tid);
			if(m==null)
				{
				m=new BinArray<T>();
				this.chroms.set(tid, m);
				}
			}
		m.put(start0, end0, o);
		return true;
		}
	
	
	public boolean isEmpty()
		{
		return this.chroms.isEmpty();
		}
	}
