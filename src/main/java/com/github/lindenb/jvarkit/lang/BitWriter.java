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
import java.io.OutputStream;

public class BitWriter {
    private final OutputStream out;
    private byte curr = '\0'; // Current byte being written
    private int offset = 0; // Tracks the position of the next bit in the byte

    public BitWriter(final OutputStream out) {
        this.out = out;
    	}

    void write(boolean bit) throws IOException {
        // Set or clear the bit at the current offset
        if (bit) {
            curr |= (1 << offset); // Set the bit
        	}
        offset++;

        // If the byte is full (8 bits written), write it to the OutputStream
        if (offset == 8) flushByte();
    	}
    
    private void flushByte() throws IOException {
    	out.write(curr);
        curr = '\0'; // Reset the current byte buffer
        offset = 0; // Reset the offset
    	}
    
    public void flush() throws IOException {
        // If there are remaining bits, flush them as a partial byte
        if (offset > 0) flushByte();
        out.flush(); // Ensure the OutputStream is flushed
    }

}
