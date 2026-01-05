/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.variant.swing;

import java.util.Vector;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import com.github.lindenb.jvarkit.swing.ColumnDef;
import com.github.lindenb.jvarkit.swing.ColumnDefTableModel;

import htsjdk.variant.variantcontext.VariantContext;

@SuppressWarnings("serial")
public class SwingVCFInfoTableModel extends ColumnDefTableModel<SwingVCFInfoTableModel.InfoLine> {
static class InfoLine {
	String key;
	Integer idx;
	Object value=null;
	}

public SwingVCFInfoTableModel() {
	super(Arrays.asList(
			new ColumnDef<InfoLine>("ID",String.class,F->F.key)	,
			new ColumnDef<InfoLine>("Index",Integer.class,F->F.idx)	,
			new ColumnDef<InfoLine>("Value",Object.class,F->F.value)	
		));
	}

public void setVariant(final VariantContext ctx) {
	setAttributes(ctx==null?null:ctx.getAttributes());
	}

public void setAttributes(final Map<String,Object> map) {
	final List<InfoLine> newrows = new Vector<>();
	if(map!=null && !map.isEmpty()) {
		for(final String key: new TreeSet<>(map.keySet()))
			{
			Object v= map.get(key);
			final List<?> L;
			if(v instanceof List)
				{
				L=(List<?>)v;
				}
			else if(v.getClass().isArray())
				{
				Object a[]=(Object[])v;
				L=Arrays.asList(a);
				}
			else
				{
				L=Collections.singletonList(v);
				}
			for(int x=0;x< L.size();++x)
				{
				final InfoLine F = new InfoLine();
				F.key = key;
				F.idx = (L.size()==1?null:x+1);
				F.value = L.get(x);
				newrows.add(F);
				}
			}
		}
	setRows(newrows);
	}


}
