/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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

import com.beust.jcommander.IStringConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;

public class FractionConverter implements IStringConverter<Double> {
	public static final String OPT_DESC="A decimal number between 0.0 and 1.0. "
			+ "If the value ends with '%' it is interpretted as a percentage eg. '1%' => '0.01'."   
			+ " A slash '/' is interpretted as a ratio. e.g: '1/100' => '0.01'."
			; 
	@Override
	public Double convert(final String s0) {
		if(StringUtils.isBlank(s0)) throw new IllegalArgumentException("Cannot convert empty string to decimal numnber.");
		String s=s0.trim();
		double value;
		try {
			if(s.endsWith("%")) {
				s=s.substring(s.length()-1).trim();
				value = Double.parseDouble(s)/100.0;
				}
			else if(s.contains("/'")) {
				final int slash = s.indexOf("/"); 
				if(slash==0 || slash+1==s.length()) throw new IllegalArgumentException("bad division in '"+s0+"'");
				final double v1 = Double.parseDouble(s.substring(0, slash));
				final double v2 = Double.parseDouble(s.substring(slash+1));
				if(v2==0) throw new IllegalArgumentException("division by zero in "+s0);
				value = v1/v2;
				}
			else
				{
				value = Double.parseDouble(s);
				}
			}
		catch(final NumberFormatException err) {
			throw new IllegalArgumentException("Cannot convert "+s0+" to decimal numnber.",err);
			}
		if(Double.isNaN(value) || Double.isInfinite(value) || value < 0.0 || value >1.0) {
			throw new IllegalArgumentException("Cannot convert "+s0+" to decimal numnber. Value should be between 0 and 1 but got '"+ value+"'.");
			}
		return value;
		}

}
