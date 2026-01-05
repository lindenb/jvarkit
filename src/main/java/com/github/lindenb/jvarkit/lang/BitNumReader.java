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

/**
 * Read integers encoded in 'n' bits
 *
 */
public class BitNumReader {
	private final double default_max;
	private final BitReader delegate;
	private final int default_nbits;
	
	private static int checkNBits(final int nbits) {
		if (nbits <= 0 || nbits > Integer.SIZE) {
		    throw new IllegalArgumentException("The number of bits (nbits) must be between 1 and "+Integer.SIZE);
			}
		return nbits;
		}
	
	private static double maxForNBits(final int nbits) {
		return Math.pow(2,nbits)-1;
		}
	
	public BitNumReader(final BitReader delegate,final int nbits) {
		this.delegate = delegate;
		this.default_nbits = checkNBits(nbits);
		this.default_max = maxForNBits(nbits);
		}
	/** get default number of bits */
	public int getDefaultNBits() {
		return this.default_nbits;
		}
	/** get default maximum value  */
	public double getDefaultMax() {
		return this.default_max;
		}
	/** get 1.0 / default maximum value  */
	public double getDefaultPrecision() {
		return 1.0 / getDefaultMax();
		}
	
	private int _nextInt(final int n_bits) throws IOException {
		int t=0;
		//int n=n_bits-1;
		for(int n=0;n< n_bits;++n) {
			switch(this.delegate.read()) {
				case 0: break;
				case 1: t+= (int)Math.pow(2,n);break;
				default: throw new IOException("cannot read "+n+"th bit");
				}
			//System.err.println("read bit["+(n_bits-n)+"] c="+n+"t0="+t+" t1="+(t*2+n));
			//MARCHE POOO t = (t << 1) | c;
			//t=t*2 + c;
			//--n;
			}
		return t;
		}
	private double _nextDouble(final int n_bits,final double max) throws IOException {
		return _nextInt(n_bits)/max;
		}
	
	
	/** get next int number [0-max] using default n_bits */
	public int nextInt() throws IOException {
		return _nextInt(getDefaultNBits());
		}
	
	/** get next number number [0-1.0] using custom 'n' bits */
	public int nextInt(int n_bits) throws IOException {
		return _nextInt(checkNBits(n_bits));
		}
	
	/** get next floating number [0-1.0] using default n_bits */
	public double nextDouble() throws IOException {
		return _nextDouble(getDefaultNBits(),getDefaultMax());
		}
	
	/** get next number number [0-1.0] using custom 'n' bits */
	public double nextDouble(int n_bits) throws IOException {
		return _nextDouble(checkNBits(n_bits),maxForNBits(n_bits));
		}
	
	@Override
	public String toString() {
		return "BitNumReader(nbits="+getDefaultNBits()+", max="+getDefaultMax()+")";
		}
	
	}
