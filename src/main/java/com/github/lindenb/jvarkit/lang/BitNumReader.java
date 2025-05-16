/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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


public class BitNumReader {
	private final double max;
	private final BitReader delegate;
	private final int nbits;
	public BitNumReader(final BitReader delegate,int nbits) {
		if (nbits <= 0 || nbits > 32) {
		    throw new IllegalArgumentException("The number of bits (nbits) must be between 1 and 32.");
			}
		this.max = Math.pow(2,nbits)-1;
		this.delegate = delegate;
		this.nbits = nbits;
		}
		
	public double next() throws IOException {
		int t=0;
		int n = this.nbits;
		while(n>0) {
			final int c = delegate.read();
			if(c==-1) throw new IOException("cannot read "+n+"th bit");
			t = (t << 1) | c;
			// t=t*2 + c;
			--n;
			}
		return t/max;
		}
	}
