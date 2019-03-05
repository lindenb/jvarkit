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
package com.github.lindenb.jvarkit.util.iterator;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;
import java.io.Reader;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.Objects;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;

/** implementation of htsjdk.tribble.readers.LineIterator */
public class LineIterator 
	extends AbstractIterator<String>
	implements htsjdk.tribble.readers.LineIterator,
	CloseableIterator<String> {

	private final Iterator<String> delegate;
	private static class BuffReadIter 
	extends AbstractIterator<String>
	implements Closeable
		{
		private BufferedReader in;
		BuffReadIter(final Reader in)
			{
			
			this.in=
					in instanceof BufferedReader?
					BufferedReader.class.cast(in):
					new BufferedReader(in);
			}
		@Override
		protected String advance()
			{
			if(in==null) return null;
			final String s;
			try {
				s = in.readLine();
				if(s==null) this.close();
				return s;
			} catch (final IOException e) {
				throw new RuntimeIOException(e);
				}
			}
		@Override
		public void close() throws IOException {
			CloserUtil.close(this.in);
			this.in=null;
			}
		}
	
	
	public LineIterator(final Iterator<String> delegate) {
		this.delegate = Objects.requireNonNull(delegate);
		}
	public LineIterator(final Collection<String> col) {
		this(Objects.requireNonNull(col).iterator());
		}
	public LineIterator(final String[] array) {
		this(Arrays.asList(array));
		}
	
	public LineIterator(final Reader br) {
		this(new BuffReadIter(br));
		}
	
	@Override
	protected String advance()
		{
		return delegate.hasNext()?delegate.next():null;
		}
	
	@Override
	public void close() {
		CloserUtil.close(this.delegate);
		}
	
	public String toString() {
		return "LineIterator";
		}
	}
