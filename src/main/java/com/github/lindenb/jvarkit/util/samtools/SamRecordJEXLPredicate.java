package com.github.lindenb.jvarkit.util.samtools;

import java.util.HashMap;
import java.util.Map;
import java.util.function.Function;
import java.util.function.Predicate;

import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlContext;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContextUtils;

public class SamRecordJEXLPredicate implements Predicate<SAMRecord> {

	private static final Map<String,Function<SAMRecord,Object>> MAPPER=new HashMap<>();
	static {
		MAPPER.put("record",R->R);
		MAPPER.put("Flags",R->R.getFlags());
		MAPPER.put("Mapped",R->R.getFlags());
		};
	
	private final String exprStr;
	private final Expression expr;
	private SamRecordJEXLPredicate(final String exprStr) {
		this.exprStr = exprStr; 
		this.expr=VariantContextUtils.engine.get().createExpression(exprStr);
	}
	
	private static class SamRecordJEXLContext
		implements JexlContext
		{
		final SAMRecord record;
	
		private SamRecordJEXLContext(final SAMRecord record) {
			this.record = record;
			}
		
		@Override
		public Object get(final String name) {
			if(MAPPER.containsKey(name)) {
				return MAPPER.get(name).apply(this.record);
				}
			return null;
			}
		@Override
		public boolean has(final String key) {
			return MAPPER.containsKey(key);
			}
		@Override
		public void set(final String key, Object arg1) {
			throw new UnsupportedOperationException();
			}
		} 
	
	
	private static final Predicate<SAMRecord> ACCEPT_ALL= new Predicate<SAMRecord>() {
		public boolean test(final SAMRecord t)
			{
			return true;
			}
		public String toString()
			{
			return "Accept all";
			}
		};
	
	public static Predicate<SAMRecord> create(final String expr) {
		if(StringUtil.isBlank(expr)) return ACCEPT_ALL;
		return new SamRecordJEXLPredicate(expr);
		}
	
	@Override
	public boolean test(final SAMRecord ctx) {
		final Object o=expr.evaluate(new SamRecordJEXLContext(ctx));
		if(o==null) return false;
		if(o instanceof Boolean) {
			return Boolean.class.cast(o).booleanValue();
		}
		if(o instanceof Integer) {
			return Integer.class.cast(o).intValue()!=0;
		}
		throw new IllegalArgumentException("expression "+this.exprStr+" doesn't return a boolean.");
		}
	
	@Override
	public String toString() {
		return this.exprStr;
		}
	
}
