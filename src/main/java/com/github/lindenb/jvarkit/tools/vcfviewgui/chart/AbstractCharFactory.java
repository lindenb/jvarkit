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
package com.github.lindenb.jvarkit.tools.vcfviewgui.chart;

import com.github.lindenb.jvarkit.tools.vcfviewgui.PedFile;


public abstract class AbstractCharFactory<HEADER,T> implements ChartFactory<HEADER, T> {
	private PedFile pedigree=PedFile.getEmptyInstance();
	private HEADER header=null;
	
	@Override
	public void setHeader(final HEADER header) {
		this.header=header;	
		}
	@Override
	public HEADER getHeader() {
		return this.header;	
		}
	public PedFile getPedigree()
		{
		return pedigree;
		}
	public void setPedigree(final PedFile pedigree)
		{
		this.pedigree = pedigree;
		}
	}
