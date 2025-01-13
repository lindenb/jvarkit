/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.bbfile;

import java.util.AbstractList;

import org.broad.igv.bbfile.BedFeature;

import htsjdk.samtools.util.Locatable;

/**
 * wraps org.broad.igv.bbfile.BedFeature as a java.util.List<String>
 *
 */
public class BigBedFeatureAsList
	extends AbstractList<String> implements Locatable {
	private final BedFeature delegate;
	private final String[] tokens;
	public BigBedFeatureAsList(final BedFeature delegate) {
		this.delegate = delegate;
		this.tokens = delegate.getRestOfFields();
		}
	@Override
	public String getContig() {return this.delegate.getChromosome();}
	@Override
	public int getStart() {return this.delegate.getStartBase()+1;}
	@Override
	public int getEnd() {return this.delegate.getEndBase();}
	@Override
	public int size() {return 3 + (this.tokens==null?0:this.tokens.length);}
	@Override
	public String get(int index)
		{
		switch(index) {
			case 0: return this.delegate.getChromosome();
			case 1: return String.valueOf(this.delegate.getStartBase());
			case 2: return String.valueOf(this.delegate.getEndBase());
			default:
				{
				index -=3;
				if(index<0 || index> tokens.length) return "";
				return tokens[index];
				}
			}
		}
	@Override
	public String toString()
		{
		return String.join("\t", this);
		}
	}
