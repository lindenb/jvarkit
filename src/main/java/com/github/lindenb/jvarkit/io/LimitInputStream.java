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
package com.github.lindenb.jvarkit.io;

import java.io.IOException;
import java.io.InputStream;

/**
 * limit input stream for security, used in VCFRecordGuesser in conjonction with pushbackinputstream
 * @author lindenb
 *
 */
public class LimitInputStream extends InputStream {
	private final InputStream delegate;
	private final long max_size;
	private long nRead = 0L;
	public LimitInputStream(final InputStream delegate,final long max_size) {
		this.delegate = delegate;
		this.max_size = max_size;
		if(this.max_size < 0L) throw new IllegalArgumentException("negative size");
		}
	@Override
	public int read() throws IOException {
		if(this.nRead>=this.max_size) return -1;
		final int c = this.delegate.read();
		if(c==-1) return -1;
		this.nRead++;
		return c;
		}
	@Override
	public int read(final byte[] b, int off, int len) throws IOException {
		long max_avail = this.max_size - this.nRead;
		if(max_avail<=0L) return -1;
		len = (int)Math.min((long)len, max_avail);
		int nRead= this.delegate.read(b, off, len);
		if(nRead==-1) return -1;
		this.nRead+= nRead;
		return nRead;
		}

	@Override
	public void close() throws IOException {
		this.delegate.close();
		this.nRead=this.max_size;
		}

	@Override
	public String toString() {
		return getClass().getName()+"("+this.max_size+")";
		}
	}
