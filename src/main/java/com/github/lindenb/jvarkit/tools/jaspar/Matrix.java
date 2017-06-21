/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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


import htsjdk.tribble.readers.LineIterator;


public class Matrix
	{
	private static final char BASES[]={'A','C','G','T'};
	public enum Type {PFM,PWM};
	private Type type;
	private String name;
	private double data[];
	
	private Matrix(String name,double data[])
		{
		this.type=Type.PFM;
		this.name=name;
		this.data=data;
		}
	
	public double score(CharSequence S)
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
		StringBuilder b=new StringBuilder(this.length());
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
		Matrix pwm=new Matrix(this.name,new double[data.length]);
		pwm.type=Type.PWM;
		
		
		/* number of observations */
		double log2= Math.log(2);
		double prior_params=0.25;
		for(int x=0;x< length();++x)
			{
			int N=(int)sum(x);
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
	
	public double get(int y,int x)
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
		StringBuilder b=new StringBuilder();
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
		implements Iterator<Matrix>
		{
		private Pattern ws=Pattern.compile("[\t ]+");
		private LineIterator iter;
		MyIterator(LineIterator iter)
			{
			this.iter=iter;
			}
		
		@Override
		public boolean hasNext()
			{
			while(iter.hasNext())
				{
				String line=iter.peek();
				if(line.isEmpty())
					{
					iter.next();
					continue;
					}
				if(!line.startsWith(">"))
					{
					throw new RuntimeException("Expected line to start with '>' but got "+line);
					}
				return true;
				}
			return false;
			}
		
		@Override
		public Matrix next() {
			if(!hasNext()) throw new IllegalStateException();
			String header=iter.next().substring(1).trim();
			double data[]=null;
			for(int i=0;i< 4;++i)
				{
				if(!iter.hasNext()) throw new IllegalStateException();
				String line=iter.next();
				String tokens[]=this.ws.split(line);
				if(i==0)
					{
					data=new double[tokens.length*4];
					}
				else
					{
					if(tokens.length*4!=data.length) throw new RuntimeException("Bad matrix in "+header);
					}
				for(int j=0;j< tokens.length;++j) data[i*tokens.length+j]=Integer.parseInt(tokens[j]);
				}
			return new Matrix(header, data);
			}
		
		@Override
		public void remove()
			{
			throw new UnsupportedOperationException();
			}
		}
	
	}
