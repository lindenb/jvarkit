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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.github.lindenb.jvarkit.swing.ColumnDef;
import com.github.lindenb.jvarkit.swing.ColumnDefTableModel;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

@SuppressWarnings("serial")
public class SwingVCFInfoHeaderLineTableModel extends ColumnDefTableModel<VCFInfoHeaderLine>{
public SwingVCFInfoHeaderLineTableModel(final VCFHeader h) {
	this(h==null?Collections.emptyList(): new ArrayList<>(h.getInfoHeaderLines()));
	}
public SwingVCFInfoHeaderLineTableModel(final List<VCFInfoHeaderLine> list) {
	super(
		Arrays.asList(
				new ColumnDef<VCFInfoHeaderLine>("ID",String.class,F->F.getID()),
				new ColumnDef<VCFInfoHeaderLine>("Type",String.class,F->F.getType()==null?null:F.getType().name()),
				new ColumnDef<VCFInfoHeaderLine>("Count",Integer.class,F->F.isFixedCount()?F.getCount():null),
				new ColumnDef<VCFInfoHeaderLine>("CountType",String.class,F->F.getCountType().name()),
				new ColumnDef<VCFInfoHeaderLine>("Description",String.class,F->F.getDescription())
				),
		list);
	}
}
