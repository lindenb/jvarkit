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
package com.github.lindenb.jvarkit.lang;

import java.io.IOException;
import java.io.InputStream;

public class BitReader {
	private final InputStream in;
	private byte curr;
	private int offset = 8;
	private final long max_to_read;
	private long nRead = 0L;
	public BitReader(final InputStream in,long max_to_read) {
		this.in = in;
		this.max_to_read=max_to_read;
		}
	public BitReader(final InputStream in) {
		this(in,-1L);
		}
	/** return 0/1 or -1 */
	int read() throws IOException {
		if(max_to_read!=-1L && nRead>=max_to_read) return -1;
		if (offset >= 8) {
			int c = in.read();
			if (c == -1)
				return -1;
			this.curr = (byte) c;
			offset = 0;
			}
		final int bit = (curr >> (offset)) & 1;
		offset++;
		nRead++;
		return bit;
		}

		
	@Override
	public String toString() {
		return "BitReader("+(int)curr+","+offset+")";
		}
	}
