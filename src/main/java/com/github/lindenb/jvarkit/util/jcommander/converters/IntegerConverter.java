package com.github.lindenb.jvarkit.util.jcommander.converters;

import java.math.BigInteger;
import java.util.function.Function;
import java.util.zip.Deflater;

import com.beust.jcommander.ParameterException;
import com.beust.jcommander.converters.BaseConverter;
import com.beust.jcommander.converters.IntegerConverter;

/**
 * More elaborate integer converter
 */
public class IntegerConverter
extends BaseConverter<Integer> implements Function<String, Integer> {
	public IntegerConverter() {
		super("");
		}
	public IntegerConverter(final String arg) {
		super(arg);
		}

	@Override
	public final Integer apply(String t) {
		return convert(t);
		}
	
	@Override
	public Integer convert(final String s) {
		final BigInteger m = new BigInteger(String.valueOf(Integer.MIN_VALUE));
		final BigInteger M = new BigInteger(String.valueOf(Integer.MAX_VALUE));
		final BigInteger n;
		
		try {
			n = new BigInteger(s);
		} catch (NumberFormatException e) {
			throw new ParameterException(getErrorString(s, "An integer"));
		}
		if(n.compareTo(m)<0)
			throw new ParameterException(getErrorString(s, "An integer >="+m));				
		if(m.compareTo(n)<0)
			throw new ParameterException(getErrorString(s, "An integer <="+M));				
		
		return n.intValueExact();
		}
	@Override
	public String toString() {
		return "Integer converter";
		}
	}
