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
package com.github.lindenb.jvarkit.iterator;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.Objects;

import com.github.lindenb.jvarkit.io.IOUtils;

import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.LineIterator;

/** implementation of htsjdk.tribble.readers.LineIterator */
public class LineIterators  {
	
	/** return LineIterators from array of Strings */
	public static LineIterator of(final String[] strings) {
		return of(Arrays.asList(strings));
		}

	/** return LineIterators from collection of Strings */
	public static LineIterator of(final Collection<String> col) {
		return of(Objects.requireNonNull(col).iterator());
		}
	
	/** return LineIterators from iterator of Strings */
	public static LineIterator of(final Iterator<String> col) {
		@SuppressWarnings("resource")
		final PeekableIterator<String> peek = new PeekableIterator<>(col);
		return new LineIterator() {
			@Override
			public boolean hasNext() {
				return peek.hasNext();
				}
			@Override
			public String next() {
				return peek.next();
				}
			@Override
			public String peek() {
				return peek.peek();
				}
			};
		}

	public static LineIterator of(final InputStream in) {
		try {
			return of(new InputStreamReader(in, "UTF-8"));
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}

	public static LineIterator of(final Reader in) {
		return of(new BufferedReader(in));
		}

	
	public static LineIterator of(final BufferedReader in) {
		return of(in.lines().iterator());
		}
	
	public static LineIterator of(final Path path) throws IOException {
		final BufferedReader br = IOUtils.openPathForBufferedReading(path);
		final LineIteratorImpl iter = new LineIteratorImpl(br.lines().iterator());
		iter.onClose = ()->{try {br.close();} catch(final Throwable err) {err.printStackTrace();}};
		return iter;
		}
	
	
	private static class LineIteratorImpl implements LineIterator {
		private final PeekableIterator<String> delegate;
		Runnable onClose = ()->{};
		LineIteratorImpl(final Iterator<String> delegate)  {
			this.delegate = new PeekableIterator<>(delegate);
			}
		@Override
		public boolean hasNext() {
			boolean b= this.delegate.hasNext();
			if(!b) onClose.run();
			return b;
			}
		@Override
		public String peek() {
			return this.delegate.peek();
			}
		@Override
		public String next() {
			return this.delegate.next();
			}
		
		}
	}
