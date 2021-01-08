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


import com.beust.jcommander.ParameterException;

/** same as RatioConverter but check it between 0 and 1 */
public class FractionConverter extends RatioConverter {
	public static final String OPT_DESC="A decimal number between 0.0 and 1.0. " + RatioConverter.OPT_DESC;
	
	@Override
	public double applyAsDouble(final String s0) {
		final double value = super.applyAsDouble(s0);
		if(Double.isNaN(value) || Double.isInfinite(value) || value < 0.0 || value >1.0) {
			throw new ParameterException("Cannot convert "+s0+" to decimal number. Value should be between 0 and 1 but got '"+ value+"'.");
			}
		return value;
		}
	}
