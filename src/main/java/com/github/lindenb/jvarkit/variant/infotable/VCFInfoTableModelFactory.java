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
package com.github.lindenb.jvarkit.variant.infotable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Collectors;

import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableModel;

import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VCFInfoTableModelFactory  {
private static final Logger LOG = Logger.of(VCFInfoTableModelFactory.class);
private static final String PREFIX="jvarkit.info.bean=";

private static class InfoTableImpl implements VCFInfoTable {
	String tag;
	String label = null;
	String desc = null;
	Function<String,List<String>> lineSplitter= S->Arrays.asList(CharSplitter.PIPE.split(S));
	final List<String> columns=new ArrayList<>();
	@Override
	public String getTag() {
		return tag;
		}
	@Override
	public List<String> getColumnNames() {
		return Collections.unmodifiableList(this.columns);
		}
	@Override
	public String getLabel() {
		return StringUtils.ifBlank(this.label, getTag());
		}
	@Override
	public String getDescription() {
		return StringUtils.ifBlank(this.desc, getLabel());
		}
	private List<String> parseOne(final String att) {
		List<String> tokens = lineSplitter.apply(att);
		// right padding with empty string
		if(tokens.size()< columns.size()) {
			tokens = new ArrayList<>(tokens);
			while(tokens.size()< columns.size())  {
				tokens.add("");
				}
			}
		return tokens;
		}

	@Override
	public List<List<String>> parse(final VariantContext ctx) {
		Objects.requireNonNull(ctx,"ctx is null");
		if(!ctx.hasAttribute(this.getTag())) return Collections.emptyList();
		final List<String> atts = ctx.getAttributeAsStringList(this.getTag(), "");
		final List<List<String>> L = new ArrayList<>(atts.size());
		for(int i=0;i< atts.size();i++) {
			L.add(parseOne(atts.get(i)));
			}
		return L;
		}
	
	@Override
	public TableModel parseAsTableModel(VariantContext ctx) {
		final List<List<String>> rows = parse(ctx);
		final DefaultTableModel tm = new DefaultTableModel(
				this.columns.stream().
					map(S->Object.class.cast(S)).
					toArray(N->new Object[N]),
				rows.size()
				);
		for(int y=0;y< rows.size();++y) {
			final List<String> row = rows.get(y);
			for(int x=0;y< row.size();++x) {
				tm.setValueAt(row.get(x), y, x);
				}
			}
		
		return tm;
		}
	@Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        final InfoTableImpl other = (InfoTableImpl) obj;
        return Objects.equals(tag, other.tag);
    	}
	@Override
	public int hashCode() {
		return tag.hashCode();
		}
	@Override
	public String toString() {
		return getTag();
		}
	}

public VCFInfoTableModelFactory() {
	
	}

/**
 * Encodes a column list as a string for the VCF header description.
 */
public static String encode(String...strings) {
	return  encode(Arrays.asList(strings));
	}


/**
 * Encodes a column list as a string for the VCF header description.
 */			

public static String encode(final List<String> columns) {
	return PREFIX+"("+String.join("|", columns)+")";
	}

/**
 * Parses a VCF INFO header line to extract a table definition.
 */
public Optional<VCFInfoTable> parseInfoHeaderLine(final VCFInfoHeaderLine infoHdrLine) {
	if(infoHdrLine==null) return Optional.empty();
	if(infoHdrLine.getType()!=VCFHeaderLineType.String) return Optional.empty();
	final String desc = infoHdrLine.getDescription();
	if(StringUtils.isBlank(desc)) return Optional.empty();
	final int i = desc.indexOf(PREFIX);
	if(i==-1 || i+PREFIX.length()==desc.length()) return Optional.empty();
	try {
		final InfoTableImpl table = new InfoTableImpl();
		table.tag = infoHdrLine.getID();
		char after = desc.charAt(i+PREFIX.length());
		if(after=='(')
			{
			
			
			final int j = desc.indexOf(')',i+PREFIX.length());
			if(j==-1) return Optional.empty();
			final List<String>  cols = table.lineSplitter.apply(desc.substring(i+PREFIX.length()+1,j).trim());
			for(int k=0;k< cols.size();k++) {
				table.columns.add(cols.get(k).trim());
				}
			
			table.desc = desc.substring(0,i)+" "+desc.substring(j+1);
			table.desc= StringUtils.normalizeSpaces(table.desc);
			}
		else if(after=='{') {
            LOG.warn("JSON-based info table encoding not yet supported for tag " + table.tag);
            
			}
		if(table.columns.isEmpty()) return Optional.empty();
		return Optional.of(table);
		}
	catch(Throwable err ) {
		return Optional.empty();
		}
	}
/**
 * Parses all INFO header lines in a VCF header, returning those that define tables.
 */
public List<VCFInfoTable> parse(final VCFHeader header) {
	return Objects.requireNonNull(header,"vcf header is null").
			getInfoHeaderLines().
			stream().
			map(H->parseInfoHeaderLine(H).orElse(null)).
			filter(H->H!=null).collect(Collectors.toList());
	}



}
