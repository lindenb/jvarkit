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
package com.github.lindenb.jvarkit.samtools.util;

import java.util.Objects;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;

public class LocatableDelegate<T extends Locatable> extends AbstractLocatable implements Comparable<LocatableDelegate<T>> {
private final T delegate;



public LocatableDelegate(final T delegate) {
	this.delegate = Objects.requireNonNull(delegate);
	}

public T getDelegate() {
	return delegate;
	}

@Override
public String getContig() {
	return getDelegate().getContig();
	}

@Override
public int getStart() {
	return getDelegate().getStart();
	}

@Override
public int getEnd() {
	return getDelegate().getEnd();
	}

@Override
public int compareTo(final LocatableDelegate<T> o) {
	return LocatableUtils.compareTo(this, o);
	}

@Override
public Interval toInterval() {
	if(getDelegate() instanceof Interval) return Interval.class.cast(getDelegate());
	return super.toInterval();
	}
}
