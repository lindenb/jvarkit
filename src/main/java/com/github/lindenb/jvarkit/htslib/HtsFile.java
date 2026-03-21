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
package com.github.lindenb.jvarkit.htslib;

import java.io.Closeable;
import java.io.Flushable;
import java.io.IOException;

import com.github.lindenb.jvarkit.jni.AbstractCPtr;


public class HtsFile extends AbstractCPtr implements Closeable, Flushable {
	public HtsFile(final String s,final String m) throws IOException {
		super(HtsLib.hts_hopen(s,m),true);
		if(isPtrNull()) throw new IOException("Cannot open "+s);
		}
	
	public boolean isOpen() {
		return isPtrNotNull();
		}
	
	
	@Override
	public void flush() throws IOException {
		if(isOpen()) HtsLib.hts_flush(getPtr());
		}
	@Override
	public final void close() {
		if(isOpen()) {
			HtsLib.hts_close(getPtr());
			setPtrToNull();
			}
		}
	@Override
	public void disposePtr() {
		if(isMemoryManaged()) close();
		}
	
	
	@Override
	public boolean equals(Object obj) {
		if(this==obj) return true;
		if(!(obj instanceof HtsFile)) return false;
		return this.getPtr()==HtsFile.class.cast(obj).getPtr();
		}
	
	
	}
