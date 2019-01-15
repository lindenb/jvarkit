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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.evs2bed;

import edu.washington.gs.evs.SnpData;
import htsjdk.tribble.Feature;

public class SnpDataFeature
implements Feature
	{
	private SnpData snpData;
	
	public SnpDataFeature(SnpData snpData)
		{
		this.snpData = snpData;
		}
	
	@Override
	@Deprecated
	public String getChr() {
		return getContig();
		}
	
	@Override
	public String getContig() {
		return getSnpData().getChromosome();
		}
	
	@Override
	public int getEnd()
		{
		return getStart();
		}
	@Override
	public int getStart() {
		return getSnpData().getChrPosition();
		}
	
	
	public SnpData getSnpData()
		{
		return snpData;
		}
	}
