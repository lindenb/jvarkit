package com.github.lindenb.jvarkit.util.jcommander.converters;

import java.util.function.Function;
import java.util.zip.Deflater;

import com.beust.jcommander.converters.IntegerConverter;

/**
 * Special converter for Zip compression. Bound the values between 0 and 9
 * "best" is interpreted as BEST_COMPRESSION
 * "none" is no compression
 */
public class CompressionConverter
extends IntegerConverter implements Function<String, Integer> {
	public CompressionConverter() {
		super("");
		}
	public CompressionConverter(final String arg) {
		super(arg);
		}

	@Override
	public final Integer apply(String t) {
		return convert(t);
		}
	
	@Override
	public Integer convert(final String s) {
		if(s!=null) {
			if(s.equals("best")) return Deflater.BEST_COMPRESSION;
			if(s.equals("none")) return Deflater.NO_COMPRESSION;
		}
		final Integer n = super.convert(s);
		if(n!=null) {
			if(n<0) return Deflater.NO_COMPRESSION;
			if(n>9) return Deflater.BEST_COMPRESSION;
		}
		return n;
	}
	@Override
	public String toString() {
		return "Compression converter";
		}
	}
