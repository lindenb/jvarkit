/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.jcommander.converter;

import java.util.function.ToDoubleFunction;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.ParameterException;

import com.github.lindenb.jvarkit.lang.StringUtils;

public class RatioConverter implements IStringConverter<Double>,ToDoubleFunction<String> {
	public static final String OPT_DESC= "If the value ends with '%' it is interpretted as a percentage eg. '1%' => '0.01'."   
			+ " A slash '/' is interpretted as a ratio. e.g: '1/100' => '0.01'."
			; 
	@Override
	public final Double convert(final String s0) {
		return this.applyAsDouble(s0);
		}
	
	@Override
	public double applyAsDouble(final String s0) {
		if(StringUtils.isBlank(s0)) throw new IllegalArgumentException("Cannot convert empty string to decimal number.");
		String s=s0.trim();
		double value;
		try {
			if(s.endsWith("%")) {
				final String s2=s.substring(0,s.length()-1).trim();
				value = Double.parseDouble(s2)/100.0;
				}
			else if(s.contains("/")) {
				final int slash = s.indexOf("/"); 
				if(slash==0 || slash+1==s.length()) throw new ParameterException("bad division in '"+s0+"'");
				final double v1 = Double.parseDouble(s.substring(0, slash));
				final double v2 = Double.parseDouble(s.substring(slash+1));
				if(v2==0) throw new ParameterException("division by zero in "+s0);
				value = v1/v2;
				}
			else
				{
				value = Double.parseDouble(s);
				}
			}
		catch(final NumberFormatException err) {
			throw new ParameterException("Cannot convert '"+s0+"' to decimal number. ("+err.getMessage()+")",err);
			}
		if(Double.isNaN(value) || Double.isInfinite(value)) {
			throw new ParameterException("Cannot convert "+s0+" to decimal number. Value is not a number or infinite '"+ value+"'.");
			}
		return value;
		}
	}
