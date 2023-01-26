/*
The MIT License (MIT)

Copyright (c) 2023 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.math;

import java.util.Objects;
import java.util.OptionalLong;
import java.util.function.LongPredicate;
import java.util.function.ToLongFunction;


public class MinMaxFilterLong implements LongPredicate {
public static final String OPT_DESC="A range of integers. Syntax (!)min-inclusive:max-inclusive. Use '!' to invert the logic: number must be outside the range. ";
private final OptionalLong minIncl;
private final OptionalLong maxIncl;
private final boolean invert;

public static class StringConverter 
	implements com.beust.jcommander.IStringConverter<MinMaxFilterLong> {
	@Override
	public MinMaxFilterLong convert(final String s) {
		return new MinMaxFilterLong(s);
		}
	}

public MinMaxFilterLong(final String s0) {
	this(s0,S->Long.parseLong(S));
	}
public MinMaxFilterLong(final String s0,final ToLongFunction<String> toLong) {
	Objects.requireNonNull(s0, "String");
	Objects.requireNonNull(toLong, "String-to-Long");
	String s=s0.trim();
	if(s.startsWith("!")) {
		this.invert = true;
		}
	else
		{
		this.invert = false;
		}
	final int colon = s.indexOf(":");
	if(colon==-1) throw new IllegalArgumentException("colon missing in "+s0);
	final String left=s.substring(0,colon).trim();
	if(left.isEmpty() || left.equals("*")) {
		this.minIncl= OptionalLong.empty();
	} else {
		this.minIncl = OptionalLong.of(toLong.applyAsLong(left));
		}
	
	final String right=s.substring(colon+1).trim();
	if(right.isEmpty() || right.equals("*")) {
		this.maxIncl= OptionalLong.empty();
	} else {
		this.maxIncl = OptionalLong.of(toLong.applyAsLong(right));
		}
	validate();
	}
public MinMaxFilterLong(final long minIncl,final long maxIncl) {
	this(minIncl,maxIncl,false);
	}
public MinMaxFilterLong(final long minIncl,final long maxIncl,boolean invert) {
	this(OptionalLong.of(minIncl),OptionalLong.of(maxIncl),invert);
	}

public MinMaxFilterLong(final OptionalLong minIncl,final OptionalLong maxIncl,boolean invert) {
	Objects.requireNonNull(minIncl, "min-including");
	Objects.requireNonNull(maxIncl, "max-including");
	this.minIncl = minIncl;
	this.maxIncl = maxIncl;
	this.invert = invert;
	validate();
	}

private void validate() {
	if(this.minIncl.isPresent() && this.maxIncl.isPresent() &&
			this.minIncl.getAsLong()> this.maxIncl.getAsLong() ) {
			throw new IllegalArgumentException(""+this.maxIncl+">"+this.maxIncl);
			}
	}

@Override
public boolean test(long i) {
	boolean keep = true;
	if(this.minIncl.isPresent() && this.minIncl.getAsLong() > i) {
		keep=false;
		}
	if(keep && this.maxIncl.isPresent() && this.maxIncl.getAsLong() < i) {
		keep=false;
		}
	return keep == !this.invert;
	}
@Override
public String toString() {
		return
			(this.invert?"!":"") + 
			(this.minIncl.isPresent()?String.valueOf(this.minIncl.getAsLong()):"*")+
			":" + 
			(this.maxIncl.isPresent()?String.valueOf(this.maxIncl.getAsLong()):"*")
			;
	}
}
