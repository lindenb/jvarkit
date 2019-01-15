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
package com.github.lindenb.jvarkit.util.vcf.swing;

import java.util.List;
import java.util.regex.Pattern;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;

/**
 * 
 * AbstractVcfTable
 *
 */
public abstract class AbstractVcfTable
	extends AbstractGenericTable<String>
	{
	private static final long serialVersionUID = 1L;
	private Pattern tab=Pattern.compile("[\t]");
	
	public abstract List<String> getSamples();
	
	public AbstractVcfTable()
		{
		}
	
	@Override
	public String getColumnName(int column)
		{
		VCFHeader.HEADER_FIELDS headers[]= VCFHeader.HEADER_FIELDS.values();
		if(column< headers.length)
			{
			return headers[column].name();
			}
		if(column==headers.length)
			{
			return "FORMAT";
			}
		column-=(headers.length+1);
		return getSamples().get(column);
		}
	
	@Override
	public int getColumnCount()
		{
		int n=VCFHeader.HEADER_FIELDS.values().length;
		if(this.getSamples().isEmpty())
			{
			return n;
			}
		return n+1+getSamples().size();
		}
	
	@Override
	public Object getValueOf(String s, int columnIndex)
		{
		if(s==null) return null;
		String tokens[]=tab.split(s,columnIndex+2);
		return columnIndex<tokens.length?tokens[columnIndex]:null;
		}
	
	public abstract AbstractVCFCodec getCodec();
	
	public VariantContext getVariantContextAt(int rowIndex)
		{
		String line = getElementAt(rowIndex);
		if(line==null) return null;
		return getCodec().decode(line);
		}
	
	
	}
