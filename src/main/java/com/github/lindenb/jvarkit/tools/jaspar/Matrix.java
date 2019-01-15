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
package com.github.lindenb.jvarkit.tools.jaspar;

import java.util.Iterator;
import java.util.regex.Pattern;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.LineIterator;

/*
example:

>MA0029.1	Mecom
A  [14 20  0 27  1 27 26  0 27  0 24 23  6 15 ]
C  [ 2  1  1  0 10  0  0  0  0  3  1  0  7  6 ]
G  [ 6  2 25  0  0  0  1 27  0  0  0  4  7  3 ]
T  [ 5  4  1  0 16  0  0  0  0 24  2  0  7  3 ]

>MA0030.1	FOXF2
A  [ 1 10 17 13  3  7  0 27 27 27  0 27 16  7 ]
C  [10  7  4  5 11  0  0  0  0  0 25  0  4  4 ]
G  [ 7  5  2  5  8 20  0  0  0  0  0  0  2  6 ]
T  [ 9  5  4  5  0  0 27  0  0  0  2  0  5 10 ]

 */
public class Matrix
	{
	private static final char BASES[]={'A','C','G','T'};
	public enum Type {PFM,PWM};
	private Type type;
	private final String name;
	private final double data[];
	
	private Matrix(final String name,double data[])
		{
		this.type=Type.PFM;
		this.name=name;
		this.data=data;
		}
	
	public double score(final CharSequence S)
		{
		if(S.length()!=this.length()) throw new IllegalArgumentException("Not same length S/matrix");
		double score=0.0;
		for (int x = 0; x < length(); x++)
			{
			int y = -1;
			switch(S.charAt(x))
				{
				case 'a':case 'A': y=0; break;
				case 'c':case 'C': y=1; break;
				case 'g':case 'G': y=2; break;
				case 't':case 'T': y=3; break;
				}
			if (y == -1) continue;
			score += get(y,x);
			}
		return score;
		}
	
	
	
	public double sum(int column)
		{
		double t=0;
		for(int y=0;y<4;++y) t+=get(y,column);
		return t;
		}
	
	public double max(int column)
		{
		double t=-1;
		for(int y=0;y<4;++y) t=Math.max(t, get(y,column));
		return t;
		}
	
	public double max()
		{
		double m=0;
		for(int x=0;x< length();++x) m+=max(x);
		return m;
		}
	
	public String getArchetype()
		{
		final StringBuilder b=new StringBuilder(this.length());
		for(int x=0;x< length();++x)
			{
			char c='A';
			double m=0;
			for(int y=0;y<4;++y)
				{
				if(m>get(y,x)) continue;
				c=BASES[y];
				m=get(y,x);
				}
			b.append(c);
			}
		return b.toString();
		}
	
	public Matrix convertToPWM()
		{
		if(type!=Type.PFM) throw new IllegalStateException();
		final Matrix pwm=new Matrix(this.name,new double[data.length]);
		pwm.type=Type.PWM;
		
		
		/* number of observations */
		double log2= Math.log(2);
		double prior_params=0.25;
		for(int x=0;x< length();++x)
			{
			final int N=(int)sum(x);
			for(int y=0;y<4;++y)
				{
				
				double f=this.get(y, x);
				double w=Math.log((f+Math.sqrt(N)*prior_params)/(N+Math.sqrt(N))/prior_params)/log2;
				pwm.data[y*length()+x]=w;
				}
			}

		return pwm;
		}
	
	public int length()
		{
		return data.length/4;
		}
	
	public double get(final int y,final int x)
		{
		return data[y*length()+x];
		}
	
	public String getName()
		{
		return name;
		}
	
	@Override
	public String toString()
		{
		final StringBuilder b=new StringBuilder();
		b.append(">").append(getName()).append('\n');
		for(int y=0;y< 4;++y)
			{
			b.append(BASES[y]);
			for(int x=0;x< length();++x)
				{
				b.append(" ").append(String.format("%3d",get(y, x)));
				}
			b.append("\n");
			}
		return b.toString();
		}

	
	/** iterator reading jaspar database */
	public static Iterator<Matrix> iterator(LineIterator r)
		{
		return new MyIterator(r);
		}
	
	private static class MyIterator
		extends AbstractIterator<Matrix>
		{
		private final Pattern ws=Pattern.compile("[\t \\[\\]]+");
		private final LineIterator iter;
		MyIterator(final LineIterator iter)
			{
			this.iter=iter;
			}
		@Override
		protected Matrix advance() {
			while(iter.hasNext())
				{
				final String L0 =iter.peek();
				if(L0.isEmpty())
					{
					iter.next();
					continue;
					}
				if(!L0.startsWith(">"))
					{
					throw new RuntimeIOException("Expected line to start with '>' but got "+L0);
					}
				final String header=iter.next().substring(1).replaceAll("[ \t]+","_").trim();
				double data[]=null;
				for(int i=0;i< 4;++i)
					{
					if(!iter.hasNext()) throw new IllegalStateException();
					final String line=iter.next();
					final String tokens[]=this.ws.split(line);
					if(i==0)
						{
						data=new double[tokens.length*4];
						}
					else
						{
						if(tokens.length*4!=data.length) throw new RuntimeException("Bad matrix in "+header);
						}
					if(!tokens[0].matches("[ATGC]")) throw new RuntimeException("line in "+line);
					for(int j=1;j< tokens.length;++j) data[i*tokens.length+j]=Integer.parseInt(tokens[j]);
					}
				return new Matrix(header, data);
				}
			return null;
			}
		}
	
	}
