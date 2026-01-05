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


public class BitNumWriter {
    private final BitWriter delegate;
    private final int nbits;
    private final double max;

    /**
     * Constructor for BitNumWriter
     * @param delegate The BitWriter instance used for writing bits
     * @param nbits The number of bits to represent each number
     */
    public BitNumWriter(BitWriter delegate, int nbits) {
        if (nbits <= 0 || nbits > 32) {
            throw new IllegalArgumentException("The number of bits (nbits) must be between 1 and 32.");
        }
        this.delegate = delegate;
        this.nbits = nbits;
        this.max = Math.pow(2, nbits) - 1; // Maximum value representable with nbits
    }

    /**
     * Writes a normalized number (between 0 and 1) as nbits.
     * @param value The number to write (must be in the range [0, 1])
     * @throws IOException If writing to the BitWriter fails
     */
    public void write(double value) throws IOException {
        if (value < 0.0 || value > 1.0) {
            throw new IllegalArgumentException("The value "+value+" must be in the range [0, 1].");
        }

        // Scale the value to an integer representation
        final int scaledValue = (int) Math.round(value * max);

        // Write the scaled value bit by bit
        for (int i = nbits - 1; i >= 0; i--) {
            boolean bit = (scaledValue & (1 << i)) != 0;
            delegate.write(bit);
        }
    }
}
