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


History:
* 2015 creation

*/
package com.github.lindenb.jvarkit.io;

import htsjdk.samtools.util.LocationAware;

import java.io.IOException;
import java.io.OutputStream;

/* implementation of LocationAware for trible indexes */
public class LocationAwareOutputStream extends OutputStream implements
		LocationAware
	{
    private final OutputStream delegate;
    private long position = 0L;

    public LocationAwareOutputStream(final OutputStream out) {
        this.delegate = out;
    }
    @Override
    public final void write(final byte[] bytes) throws IOException {
        write(bytes, 0, bytes.length);
    }
    @Override
    public final void write(final byte[] bytes, final int startIndex, final int numBytes) throws IOException {
        position += numBytes;
        delegate.write(bytes, startIndex, numBytes);
    }
    @Override
    public final void write(final int c)  throws IOException {
        position++;
        delegate.write(c);
    }

    @Override
    public final long getPosition() { return position; }

    @Override
    public void close() throws IOException {
        super.close();
        delegate.close();
    }


}
