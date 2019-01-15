/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.util.align;

import java.util.Arrays;

import com.github.lindenb.jvarkit.lang.SubSequence;


public class LongestCommonSequence 
	{
	public static class Hit
		{
		private int score=0;
		private int x=0;
		private int y=0;
		private CharSequence a;
		private CharSequence b;
		
		public CharSequence getMatchingSequence()
			{
			return new SubSequence(getX(), getStartX(), getEndX());
			}
		
		@Override
		public String toString()
			{
			StringBuilder b=new StringBuilder(size());
			for(int i=0; i < size() ; ++i)
	    		{
	    		char base= getX().charAt(getStartX()+i);
	    		b.append(base);
	    		}
			return b.toString();
			}
		public int size()
			{
			return this.score;
			}
		
		public CharSequence getX()
			{
			return a;
			}
		
		public CharSequence getY()
			{
			return b;
			}
		
		public int getStartX()
			{
			return x;
			}
		public int getStartY()
			{
			return y;
			}
		public int getEndX()
			{
			return getStartX()+size();
			}
		public int getEndY()
			{
			return getStartY()+size();
			}
		}

	private int width;
	private int height;
	private int[] array=null;

	public boolean compare(char c1,char c2)
		{
		return c1!='N' && c1==c2;
		}
	
	public Hit align(
			    final CharSequence S1,
			    int start1,
			    int end1,
				final CharSequence S2,
				int start2,
				int end2
				)
		{
		final int L1=(end1-start1);
		final int L2=(end2-start2);
		/** resize matrix */
		this.width = L1+1;
		this.height = L2+1;
		if( array==null || array.length< (width*height) )
			{
			this.array = new int[width*height];
			}
		/** reset matrix */
		Arrays.fill(this.array, 0);
		
		int best_x=0;
		int best_y=0;
		int max_score=0;
		this.array[0]=0;
		for(int x=0;x< L1 ;++x)
			{
			this.array[x+1]=0;
			}
		for(int y=0;y< L2 ;++y)
			{
			this.array[(y+1)*width]=0;
			for(int x=0;x< L1 ;++x)
				{
				char c1 = S1.charAt(start1 + x);
				char c2 = S2.charAt(start2 + y);
				int v;
				if( compare(c1,c2) )
					{
					v = 1 + 
						this.array[(y)*width+(x)]//diagonal
						;
					}
				else
					{
					v = 0;
					}
				this.array[(y+1)*width+(x+1)] = v;
				if(v>max_score)
					{
					best_x  = x;
					best_y  = y;
					max_score=v;
					}
				}

			}
		Hit hit= new Hit();
		hit.a  = S1;
		hit.b  = S2;
		hit.score = max_score;
		hit.x = start1 + best_x - (max_score-1);
		hit.y = start2 + best_y - (max_score-1);
		return hit;
		}
	
}
