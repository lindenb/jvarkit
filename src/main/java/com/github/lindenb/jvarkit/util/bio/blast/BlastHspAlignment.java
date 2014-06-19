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
			return "align["+getIndex()+"]="+getQueryChar()+getMidChar()+getHitChar()+
					"/"+
					this.getQueryIndex1()+","+this.getHitIndex1();
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
	
	
	public static Hsp cloneHsp(Hsp hsp)
		{
		Hsp h2=new Hsp();
		h2.setHspNum(hsp.getHspNum());
		h2.setHspBitScore(hsp.getHspBitScore());
		h2.setHspScore(hsp.getHspScore());
		h2.setHspEvalue(hsp.getHspEvalue());
		h2.setHspIdentity(hsp.getHspIdentity());
		h2.setHspGaps(hsp.getHspGaps());
		h2.setHspPositive(hsp.getHspPositive());
		h2.setHspAlignLen(hsp.getHspAlignLen());
		
		
		h2.setHspQueryFrom(hsp.getHspQueryFrom());
		h2.setHspQueryTo(hsp.getHspQueryTo());
		h2.setHspQueryFrame(hsp.getHspQueryFrame());
		
		
		h2.setHspHitFrom(hsp.getHspHitFrom());
		h2.setHspHitTo(hsp.getHspHitTo());
		h2.setHspHitFrame(hsp.getHspHitFrame());

		h2.setHspQseq(hsp.getHspQseq());
		h2.setHspHseq(hsp.getHspHseq());
		h2.setHspMidline(hsp.getHspMidline());
		
		h2.setHspPatternFrom(hsp.getHspPatternFrom());
		h2.setHspPatternTo(hsp.getHspPatternTo());
		
		return h2;
		}
	
	
	private int alignLen=-1;
	public int getAlignLength()
		{
		if(this.alignLen<0)
			{
			if(getHsp().getHspAlignLen()!=null)
				{
				this.alignLen=Integer.parseInt(getHsp().getHspAlignLen());
				}
			else
				{
				this.alignLen=getHsp().getHspMidline().length();
				}
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
			curr.indexInAlignment=0;
			}
		
		BlastHspAlignment owner()
			{
			return BlastHspAlignment.this;
			}
		
		@Override
		public boolean hasNext()
			{
			return curr.indexInAlignment< owner().getAlignLength();
			}
		@Override
		public Align next()
			{
			if(curr.indexInAlignment >= owner().getAlignLength())
				{
				throw new IllegalStateException();
				}
			
			Align ret= curr.clone();
			
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

			
			++curr.indexInAlignment;
			return ret;
			}
		
		@Override
		public void remove()
			{
			throw new UnsupportedOperationException();
			}
		}
	}
