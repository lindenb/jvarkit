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


import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public abstract class VariantContextChartFactory
extends AbstractCharFactory<VCFHeader,VariantContext>{

	protected List<Double> getAttributeAsDoubleList(final VariantContext ctx,final String att)
		{
		if(ctx==null || att==null) return Collections.emptyList();	
		final List<Object> L=ctx.getAttributeAsList(att);
		final List<Double> F=new ArrayList<>(L.size());
		for(final Object o: L) 
			{
			if(o==null) continue;
			final double v;
			if(o instanceof Double)
				{
				v=Double.class.cast(o);
				}
			else if(o instanceof Float)
				{
				v=Float.class.cast(o);
				}
			else
				{
				try
					{
					v=Double.parseDouble(o.toString());
					}
				catch(NumberFormatException err) {
					continue;
					}
				}
			F.add(v);
			}
		return F;
		}

}
