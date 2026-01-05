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
package com.github.lindenb.jvarkit.io;

import java.io.IOException;
import java.io.Writer;
import java.util.Objects;
import java.util.function.Supplier;

import htsjdk.samtools.util.RuntimeIOException;

/**
 * A writing adding suffix of prefix after/before each carriage return.
 *
 */
public class PrefixSuffixWriter extends Writer {
	private Writer delegate;
	private boolean at_start = true;
	private Supplier<String> prefixSupplier = ()->null;
	private Supplier<String> suffixSupplier = ()->null;
	
	public PrefixSuffixWriter(Writer delegate) {
		this.delegate = Objects.requireNonNull(delegate);
	}
	
	public PrefixSuffixWriter(final String prefix,Writer delegate,String suffix) {
		this.delegate = Objects.requireNonNull(delegate);
		setPrefix(prefix);
		setSuffix(prefix);
	}
	
	public PrefixSuffixWriter(final Supplier<String> prefix,Writer delegate,Supplier<String> suffix) {
		this.delegate = Objects.requireNonNull(delegate);
		setPrefix(prefix);
		setSuffix(prefix);
	}
	
	public PrefixSuffixWriter setPrefix(Supplier<String> supplier) {
		this.prefixSupplier = (supplier==null?()->null:supplier);;
		return this;
		}
	public PrefixSuffixWriter setPrefix(final String prefix) {
		return setPrefix(()->prefix);
		}
	
	public PrefixSuffixWriter setSuffix(Supplier<String> supplier) {
		this.suffixSupplier = (supplier==null?()->null:supplier);
		return this;
		}
	public PrefixSuffixWriter setSuffix(final String prefix) {
		return setSuffix(()->prefix);
		}
	
	public PrefixSuffixWriter setNoPrefix() {
		this.prefixSupplier = ()->null;
		return this;
		}
	public PrefixSuffixWriter setNoSuffix() {
		this.suffixSupplier = ()->null;
		return this;
		}
	
	public final void println(final String s)  {
		print(s);
		println();
		}
	public final void print(final String s)  {
		try {
			write(s);
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	public final void println()  {
		print("\n");
		}

	@Override
	public void write(char[] cbuf, int off, int len) throws IOException {
		final char[] buffer1 = new char[1];
		for(int i=0;i< len;i++) {
			if(at_start) {
				final String prefix = prefixSupplier.get();
				if(prefix!=null && !prefix.isEmpty()) this.delegate.write(prefix);
				at_start=false;
				}
			final char c= cbuf[off+i];
			if(c=='\n') {
				final String suffix = suffixSupplier.get();
				if(suffix!=null && !suffix.isEmpty()) this.delegate.write(suffix);
				at_start = true;
				}
			buffer1[0]=c;
			delegate.write(buffer1, 0, 1);
			}
		}

	@Override
	public void flush() throws IOException {
		this.delegate.flush();
		
	}

	@Override
	public void close() throws IOException {
		this.delegate.close();
	}

}
