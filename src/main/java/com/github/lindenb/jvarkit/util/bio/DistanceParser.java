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
package com.github.lindenb.jvarkit.util.bio;

import java.math.BigInteger;
import java.util.function.ToIntFunction;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.ParameterException;


public class DistanceParser implements ToIntFunction<String> {
	public static final String OPT_DESCRIPTION =
			"A distance specified as a positive integer."
			+ "Comma are removed. "
			+ "The following suffixes are interpreted : b,bp,k,kb,m,mb";

	
	public static class StringConverter
		implements IStringConverter<Integer>
		{
		@Override
		public Integer convert(final String distStr) {
			try {
				final DistanceParser parser = new DistanceParser();
			      return parser.applyAsInt(distStr);
			    } catch(final Exception ex) {
			      throw new ParameterException("Bad distance \""+distStr+"\"",ex);
			    }
			}
		}
	
	private BigInteger _parseBigInteger( String s)
		{
		BigInteger factor = BigInteger.ONE;
		s=s.replace(",", "");
		
		 if(s.endsWith("bp"))
			{
			s=s.substring(0, s.length()-2).trim();
			}
		else if(s.toLowerCase().endsWith("kb"))
			{
			s=s.substring(0, s.length()-2).trim();
			factor = BigInteger.valueOf(1_000L);
			}
		else if(s.toLowerCase().endsWith("mb"))
			{
			s=s.substring(0, s.length()-2).trim();
			factor = BigInteger.valueOf(1_000_000L);
			}
		else if(s.endsWith("b"))
			{
			s=s.substring(0, s.length()-1).trim();
			}
		else if(s.endsWith("k"))
			{
			s=s.substring(0, s.length()-1).trim();
			factor = BigInteger.valueOf(1_000L);
			}
		else if(s.endsWith("m"))
			{
			s=s.substring(0, s.length()-1).trim();
			factor = BigInteger.valueOf(1_000_000L);
			}
		final BigInteger vbi = new BigInteger(String.valueOf(s)).multiply(factor);
		
		return vbi;
		}
	
	/** parse distance as BigInteger */
	public BigInteger parseAsBigInteger(final String distance) {
		try {
			return _parseBigInteger(distance);
			}
		catch(final Throwable err)
			{
			throw new IllegalArgumentException("Cannot parse distance \""+distance+"\"",err);
			}
		}
	@Override
	public int applyAsInt(final String distance) {
		final BigInteger bi = _parseBigInteger(distance);
		if(bi.compareTo(BigInteger.ZERO)<0) {
			throw new IllegalArgumentException("Cannot parse negative distance \""+distance+"\". It must be greater than zero");
			}
		try {
			return bi.intValueExact();
			}
		catch(final ArithmeticException err) {
			throw new IllegalArgumentException("Cannot parse distance \""+distance+"\"",err);
			}
		}
}
