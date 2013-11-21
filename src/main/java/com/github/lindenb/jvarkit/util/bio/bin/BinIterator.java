package com.github.lindenb.jvarkit.util.bio.bin;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/** memory safe ucsc-bin iterator */
public class BinIterator implements Iterator<Integer>
	{
	private final int MAX_GENOME_SIZE=536870912;//2^29
	private boolean _has_next_called=false;
	private boolean _has_next=false;
	private int _next=-1;
	private int beg;
	private int end;
	private int depth=0;
	private int k=0;
	public BinIterator(int beg0,int end0)
		{
		if (end0 >= MAX_GENOME_SIZE) end0 = MAX_GENOME_SIZE;
		--end0;
		this.beg=beg0;
		this.end=end0;
		}		
	
	private boolean shift()
		{
		if(beg>end) return false;
		for(;;)
			{
			switch(depth)
				{
				case 0:
					_next=0;
					depth++;
					k= 1 + (beg>>26);
					return true;
				case 1:
					if(k<=1 + (end>>26))
						{
						_next=k;
						++k;
						return true;
						}
					else
						{
						k=9 + (beg>>23);
						++depth;
						}
					break;
				case 2:
					if(k<=  9 + (end>>23))
						{
						_next=k;
						++k;
						return true;
						}
					else
						{
						k= 73 + (beg>>20);
						++depth;
						}
					break;
				case 3:
					if(k<= 73 + (end>>20))
						{
						_next=k;
						++k;
						return true;
						}
					else
						{
						k= 585 + (beg>>17);
						++depth;
						}
					break;
				case 4:
					if(k<= 585 + (end>>17))
						{
						_next=k;
						++k;
						return true;
						}
					else
						{
						k= 4681 + (beg>>14);
						++depth;
						}
					break;
				case 5:
					if(k<= 4681 + (end>>14))
						{
						_next=k;
						++k;
						return true;
						}
					else
						{
						++depth;
						return false;
						}
				default: return false;
				}
			}
		}
		
 	public boolean hasNext()
	 	{
		if(_has_next_called) return _has_next;
		_has_next_called=true;
		_has_next=this.shift();
		return _has_next;
	 	}
	public Integer next()
		{
		if(!_has_next_called) hasNext();
		if(!_has_next) throw new IllegalStateException();
		int tmp=_next;
		_next=-1;
		_has_next_called=false;
		_has_next=false;
		return tmp;
		}
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
		}
	
	@Override
	public String toString() {
		return "bin("+beg+"-"+(end+1)+")";
		}
	
	static public List<Integer> bins(int beg,int end)
		{
		List<Integer> L=new ArrayList<Integer>();
		BinIterator iter=new BinIterator(beg, end);
		while(iter.hasNext())
			{
			L.add(iter.next());
			}
		return L;
		}
	public static void main(String[] args)
		{
		if(args.length!=2)
			{
			System.err.println("Usage: begin end");
			System.exit(-1);
			}
		try
			{
			BinIterator iter=new BinIterator(
					Integer.parseInt(args[0]),
					Integer.parseInt(args[1]));
			while(iter.hasNext())
				{
				System.out.println(iter.next());
				}
			}
		catch(Exception err)
			{
			err.printStackTrace();
			System.exit(-1);
			}
		}
 	}
