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
package com.github.lindenb.jvarkit.util.samtools;

import java.util.HashMap;
import java.util.Map;
import java.util.function.Function;

import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlContext;
import org.apache.commons.jexl2.JexlException;

import com.beust.jcommander.IStringConverter;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContextUtils;

/**
 * 
 * SamRecordJEXLFilter
 *
 */
public class SamRecordJEXLFilter
	implements SamRecordFilter {
	public static final String FILTER_DESCRIPTION = 
			"A JEXL Expression that will be used to filter out some sam-records (see https://software.broadinstitute.org/gatk/documentation/article.php?id=1255). "+
			"An expression should return a boolean value (true=exclude, false=keep the read). "
			+ "An empty expression keeps everything. "
			+ "The variable 'record' is the current observed read, an instance of SAMRecord (https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html)."
			;
    public static final String DEFAULT_FILTER = 
    		  "record.getMappingQuality()<1 || "
    		+ "record.getDuplicateReadFlag() || "
    		+ "record.getReadFailsVendorQualityCheckFlag() || "
    		+ "record.isSecondaryOrSupplementary()";
    public static final String DEFAULT_OPT = "--samFilter";

	private static final Map<String,Function<SAMRecord,Object>> MAPPER=new HashMap<>();
	static {
		MAPPER.put("record",R->R);
		MAPPER.put("CLIPPED",R->{
			if(R.getReadUnmappedFlag()) return false;
			final Cigar c=R.getCigar();
			return(c!=null && c.isClipped());
		});
		};
	
	private final String exprStr;
	private final Expression expr;
	
	private SamRecordJEXLFilter(final String exprStr) {
		this.exprStr = exprStr; 
		try {
			this.expr=VariantContextUtils.engine.get().createExpression(exprStr);
		} catch(final JexlException err) {
			throw new IllegalArgumentException("Cannot compile JEXL expression", err);
		}
	}
	
	
	public static class StringConverter
	implements IStringConverter<SamRecordFilter>
		{
		@Override
		public SamRecordFilter convert(final String expr) {
			return create(expr);
			}	
		}

	private static class SamRecordJEXLContext
		implements JexlContext
		{
		final SAMRecord record;

		private SamRecordJEXLContext(final SAMRecord record) {
			this.record = record;
			if(this.record==null) throw new NullPointerException("SAMRecord is null");
			}
		@Override
		public Object get(final String name) {
			final Function<SAMRecord,Object> fun = MAPPER.get(name);
			return (fun!=null?fun.apply(this.record):null);
			}
		@Override
		public boolean has(final String key) {
			return MAPPER.containsKey(key);
			}
		@Override
		public void set(final String key, Object arg1) {
			throw new UnsupportedOperationException();
			}
		@Override
		public String toString() {
			return "JexlContext for SAMRecord "+this.record;
			}
		} 
	
	
	private static final SamRecordFilter ACCEPT_ALL= new SamRecordFilter() {
		@Override
		public boolean filterOut(final SAMRecord record) {
			return false;
			}
		@Override
		public boolean filterOut(final SAMRecord first,final SAMRecord second) {
			return false;
			}
		public String toString()
			{
			return "'Accept all' (Empty expression)";
			}
		};
	
   public static SamRecordFilter buildDefault() {
            return create(DEFAULT_FILTER);
            }

    public static SamRecordFilter buildAcceptAll() {
            return ACCEPT_ALL;
            }

		
	public static SamRecordFilter create(final String expr) {
		if(StringUtil.isBlank(expr)) return buildAcceptAll();
		return new SamRecordJEXLFilter(expr);
		}
	
	@Override
	public boolean filterOut(final SAMRecord record) {
		final Object o;
		try {
			o = this.expr.evaluate(new SamRecordJEXLContext(record));
			}
		catch(final JexlException err) {
			throw new RuntimeException("Cannot evaluate JEXL expression \""+this.exprStr+"\" with SAMRecord 'record' :"+record, err);
			}
		
		if(o==null) return true;
		if(o instanceof Boolean) {
			return Boolean.class.cast(o).booleanValue();
		}
		if(o instanceof Integer) {
			return Integer.class.cast(o).intValue()!=0;
		}
		throw new IllegalArgumentException("expression "+this.exprStr+" doesn't return a boolean.");
		}
	
	@Override
	public boolean filterOut(final SAMRecord first,final SAMRecord second) {
		return filterOut(first) && filterOut(second);
		}
	
	
	@Override
	public String toString() {
		return this.exprStr;
		}
	
}
