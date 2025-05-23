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
package com.github.lindenb.jvarkit.variant.infotable;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import javax.swing.table.TableModel;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * manage data in VCF INFO column with separated data like "VEP" "SNPEFF" etc...
 * used to display data like a table
 *
 */
public interface VCFInfoTable  {
	public String getTag();
	public String getLabel();
	public String getDescription();
	public List<String> getColumnNames();
	 /**
     * Parses the INFO field values into a list of rows (each a list of columns).
     */
	public List<List<String>> parse(final VariantContext ctx);
	/**
     * Converts INFO field values into a Swing TableModel.
     */
	public TableModel parseAsTableModel(final VariantContext ctx);
	/**
     * Parses the INFO field values into a Map of rows (each a list of columns).
     */
	public default List<Map<String,String>> parseAsMap(final VariantContext ctx) {
		final List<List<String>> L= parse(ctx);
		final List<Map<String,String>> retList = new ArrayList<>(L.size());
		final List<String> colNames = getColumnNames();
		for(int x=0;x<L.size();++x) {
			final List<String> row = L.get(x);
			final Map<String,String> map = new LinkedHashMap<>(colNames.size());
			for(int i=0;i< colNames.size();i++) {
				map.put(colNames.get(i), i< row.size()?row.get(i):"");
				}
			retList.add(map);
			}
		return retList;
	}
}
