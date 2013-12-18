package com.github.lindenb.jvarkit.util.bio.blast;

import java.util.Iterator;

import gov.nih.nlm.ncbi.blast.Hsp;

public class BlastHspAlignment
	implements Iterable<BlastHspAlignment.Align>
	{
	private final Hsp _hsp;
	
	public enum CigarOperator
		{
		EQ,X,I,D;
		}
	
	public class Align
		implements Cloneable
		{
		private int indexInAlignment=0;
		private int query_index ;
		private int hit_index;

		private Align()
			{
			}
		
		public int getQueryIndex1()
			{
			return query_index;
			}
		public int getHitIndex1()
			{
			return hit_index;
			}

		
		public Hsp getHsp()
			{
			return BlastHspAlignment.this.getHsp();
			}
		
		public int getIndex()
			{
			return this.indexInAlignment;
			}
		public char getHitChar()
			{
			return getHsp().getHspHseq().charAt(getIndex());
			}
		
		public char getQueryChar()
			{
			return getHsp().getHspQseq().charAt(getIndex());
			}

		public char getMidChar()
			{
			return getHsp().getHspMidline().charAt(getIndex());
			}
		
		public CigarOperator getCigarOperator()
			{
			//deletion in REF, insert in READ
			if(getHitChar()=='-')
				{
				return CigarOperator.I;
				}
			//deletion in READ, insert in REF
			else if(getQueryChar()=='-')
				{
				return CigarOperator.D;
				}
			else if(getMidChar()=='|')
				{
				return CigarOperator.EQ;
				}
			else
				{
				return CigarOperator.X;
				}
			}
		
		@Override
		protected Align clone()
			{
			Align a=new Align();
			a.indexInAlignment=indexInAlignment;
			a.hit_index=hit_index;
			a.query_index=query_index;
			return a;
			}
		@Override
		public boolean equals(Object obj)
			{
			if(obj==this) return true;
			if(obj==null || !(obj instanceof Hsp)) return false;
			Align cp=Align.class.cast(obj);
			return cp.getHsp()==getHsp() && cp.getIndex()==getIndex();
			}
		@Override
		public int hashCode()
			{
			return indexInAlignment;
			}
		
		@Override
		public String toString()
			{
			return "align["+getIndex()+"]";
			}
		}
	public BlastHspAlignment(final Hsp hsp)
		{
		if(hsp==null) throw new NullPointerException();
		this._hsp=hsp;
		}
	
	public Hsp getHsp()
		{
		return _hsp;
		}
	
	private int alignLen=-1;
	public int getAlignLength()
		{
		if(this.alignLen<0)
			{
			this.alignLen=Integer.parseInt(getHsp().getHspAlignLen());
			}
		return this.alignLen; 
		}
	
	private int queryFrom=-1;
	public int getQueryFrom1()
		{
		if(this.queryFrom<0)
			{
			this.queryFrom=Integer.parseInt(getHsp().getHspQueryFrom());
			}
		return this.queryFrom; 
		}
	
	private int queryTo=-1;
	public int getQueryTo1()
		{
		if(this.queryTo<0)
			{
			this.queryTo=Integer.parseInt(getHsp().getHspQueryTo());
			}
		return this.queryTo; 
		}
	
	private int hitFrom=-1;
	public int getHitFrom1()
		{
		if(this.hitFrom<0)
			{
			this.hitFrom=Integer.parseInt(getHsp().getHspHitFrom());
			}
		return this.hitFrom; 
		}
	
	private int hitTo=-1;
	public int getHitTo1()
		{
		if(this.hitTo<0)
			{
			this.hitTo=Integer.parseInt(getHsp().getHspHitTo());
			}
		return this.hitTo; 
		}

	
	private int hitShift()
		{
		return getHitFrom1()>getHitTo1()?-1:1;
		}
	
	public boolean isPlusPlus()
		{
		return hitShift()==1;
		}
	
	public final boolean isPlusMinus()
		{
		return !isPlusPlus();
		}
	
	@Override
	public Iterator<Align> iterator()
		{
		return new MyIterator();
		}
	private class MyIterator
		implements Iterator<BlastHspAlignment.Align>
		{
		private Align curr;
		MyIterator()
			{
			curr=new Align();
			curr.query_index = owner().getQueryFrom1();
			curr.hit_index = owner().getHitFrom1();
			}
		
		BlastHspAlignment owner()
			{
			return BlastHspAlignment.this;
			}
		
		@Override
		public boolean hasNext()
			{
			return curr.indexInAlignment+1< owner().getAlignLength();
			}
		@Override
		public Align next()
			{
			if(curr.indexInAlignment+1 >= owner().getAlignLength())
				{
				throw new IllegalStateException();
				}
			++curr.indexInAlignment;
			
			char ch= curr.getHitChar();
			char cq= curr.getQueryChar();

			if(ch!='-' && ch!=' ')
				{
				curr.hit_index+=hitShift();
				}
			if(cq!='-' && cq!=' ')
				{
				curr.query_index++;
				}

			return curr.clone();
			}
		
		@Override
		public void remove()
			{
			throw new UnsupportedOperationException();
			}
		}
	}
