package com.github.lindenb.jvarkit.math;

import java.util.Arrays;
import java.util.Iterator;


/* see also http://trove4j.sourceforge.net/javadocs/gnu/trove/list/array/TDoubleArrayList.html */
public class DoubleArray 
	implements Iterable<Double>,Cloneable
	{
	private double _array[];
	private int _len=0;
	private int _extend=-1;
	
	public DoubleArray(final DoubleArray copy)
		{
		this._len=copy._len;
		this._array=Arrays.copyOf(copy._array, copy._len);
		this._extend=copy._extend;
		}
	public DoubleArray(int capacity)
		{
		this(capacity,-1);
		}
	public DoubleArray(int capacity,int extend)
		{
		this._array=new double[capacity];
		this._extend=extend;
		Arrays.fill(this._array, 0);
		}
	
	public void push_back(double v)
		{
		if(_len>=this._array.length)
			{
			int x=this._extend;
			if(x==-1) x=1+(_array.length)/2;
			this._array=Arrays.copyOf(_array, _len+x);
			}	
		this._array[_len]=v;
		++_len;
		}
	
	public boolean isEmpty()
		{
		return _len==0;
		}
	
	public int size()
		{
		return _len;
		}
	
	public double getQuick(int index)
		{
		return _array[index];
		}
	
	public void setSize(int newLen,int fillValue)
		{
		if(newLen< this._len)
			{
			_len=newLen;
			}
		else if(newLen> this._len)
			{
			while(size()< newLen) push_back(fillValue);
			}
		}
	
	private void ensureIndex(int index)
		{
		if(index<0 || index >=_len) throw new IndexOutOfBoundsException("0<="+index+"<"+_len);
		}
	
	public double get(int index)
		{
		ensureIndex(index);
		return getQuick(index);
		}
	
	public double setQuick(int index,double v)
		{
		double old= _array[index];
		_array[index]=v;
		return old;
		}
	
	public double set(int index,double v)
		{
		ensureIndex(index);
		return setQuick(index,v);
		}
	@Override
	public Iterator<Double> iterator()
		{
		return new MyIterator();
		}
	@Override
	public DoubleArray clone()
		{
		return new DoubleArray(this);
		}
		
	
	public double getMin()
		{
		double minDepth=Double.MAX_VALUE;
		for(int i=0;i< _len;++i)
			{
			double v=_array[i];
			if(Double.isNaN(v)) continue;
			minDepth=Math.min(minDepth,v);
			}
		return minDepth;
		}
	
	public void substract(double v)
		{
		for(int i=0;i< this._len;++i)
			{
			this._array[i]-=v;
			}
		}
	
	public void divide(double v)
		{
		if(v==0) throw new IllegalArgumentException();
		for(int i=0;i< this._len;++i)
			{
			this._array[i]/=v;
			}
		}
	
	public void sort()
		{
		Arrays.sort(this._array,0,this._len);
		}
	
	public boolean isSorted()
		{
		if(_len<=1) return true;
		for(int i=1;i< _len;++i)
			{
			if(_array[i-1]>_array[i]) return false;
			}
		return true;
		}
	
	  public double[] toArray()
	  	{
		  double cp[]=new double[size()];
		  System.arraycopy(this._array, 0, cp, 0, this._len);
		  return cp;
	  	}
	
	
	private class MyIterator implements Iterator<Double>
		{
		private int index=0;
		@Override
		public boolean hasNext() {
			return index < DoubleArray.this.size();
			}
		@Override
		public Double next()
			{
			return DoubleArray.this.get(index++);
			}
		
		@Override
		public void remove()
			{
			throw new UnsupportedOperationException();
			}
		}
	
	  /** {@inheritDoc} */
    @Override
    public boolean equals( Object o )
    	{
    	if(o==null || !(o instanceof DoubleArray)) return false;
        if ( o == this )  return true;
        
        DoubleArray other=DoubleArray.class.cast(o);
        if ( other.size() != this.size() ) return false;
        for ( int i = 0; i < _len; i++ )
              {
               if ( this._array[ i ] != other._array[ i ] ) return false;
               }
        return true;
        }


    /** {@inheritDoc} */
    @Override
    public int hashCode()
    	{
        int h = 0;
        for ( int i = 0; i < _len; i++ )
        	{
            h += 31*((Double) this._array[ i ]).hashCode();
        	}
        return h;
    	}


    // procedures
	
	 public String toString()
 		{
        final StringBuilder buf = new StringBuilder( "[" );
        for ( int i = 0; i < _len; i++ )
        	{
        	if(i>0) buf.append( "," );
            buf.append( _array[ i ] );
        	}
        buf.append( "]" );
        return buf.toString();
	    }


	}
