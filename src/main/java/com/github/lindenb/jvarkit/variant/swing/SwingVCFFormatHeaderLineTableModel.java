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

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.swing.ColumnDef;
import com.github.lindenb.jvarkit.swing.ColumnDefTableModel;

import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

@SuppressWarnings("serial")
public class SwingVCFFormatHeaderLineTableModel extends ColumnDefTableModel<VCFFormatHeaderLine>{
public SwingVCFFormatHeaderLineTableModel(final List<VCFFormatHeaderLine> list) {
	super(Arrays.asList(
		 new ColumnDef<VCFFormatHeaderLine>("ID",String.class,F->F.getID()),
		 new ColumnDef<VCFFormatHeaderLine>("Type",String.class,F->F.getType()==null?null:F.getType().name()),
		 new ColumnDef<VCFFormatHeaderLine>("Count",Integer.class,F->F.isFixedCount()?F.getCount():null),
		 new ColumnDef<VCFFormatHeaderLine>("CountType",String.class,F->F.getCountType().name()),
		 new ColumnDef<VCFFormatHeaderLine>("Description",String.class,F->F.getDescription())
		),list);
	}

public SwingVCFFormatHeaderLineTableModel(final VCFHeader header) {
	this(header.getFormatHeaderLines().stream().
			sorted((A,B)->A.getID().compareTo(B.getID())).
			collect(Collectors.toList()));
	}

}
