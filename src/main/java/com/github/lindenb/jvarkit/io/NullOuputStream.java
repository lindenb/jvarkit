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
package com.github.lindenb.jvarkit.io;

import java.io.IOException;
import java.io.OutputStream;
/** output stream that doesn't print anything */
public class NullOuputStream extends OutputStream
	{
	private boolean _closed=false;
	private long count=0L;
	public NullOuputStream()
		{
		}
	@Override
	public void write(int c)// throws IOException
		{
		if(!isClosed() && c>=0) count++;
		}
	@Override
	public void write(byte[] b)// throws IOException
		{
		if(!isClosed()) count+=b.length;
		}
	@Override
	public void write(byte[] b, int off, int len)// throws IOException
		{
		if(!isClosed()) count+=len;
		}
	@Override
	public void close() throws IOException
		{
		_closed=true;
		}
	public long getByteWrittenCount()
		{
		return count;
		}
	public boolean isClosed()
		{
		return _closed;
		}
	
	@Override
	public void flush()  {
		}
	@Override
	public String toString()
		{
		return "NullOuputStream: closed:"+ isClosed()+" N="+getByteWrittenCount();
		}
	}
