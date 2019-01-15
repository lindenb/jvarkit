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
package com.github.lindenb.jvarkit.htslib;

import java.io.Closeable;
import java.io.IOException;

import com.github.lindenb.jvarkit.jni.CPtr;

public class HtsFile extends CPtr implements Closeable {
	public HtsFile(final String s,final String m) throws IOException {
		super(_open(s,m));
		if(isNull()) throw new IOException("Cannot open "+s);
		}
	@Override
	public void close() {
		if(!isNull()) _close(this.get());
		setNull();
		}
	
	public boolean readLine(char delim,final KString ks)  throws IOException {
		return _getline(get(),delim,ks.get())!=-1;
		}
	
	@Override
	public void dispose() {
		close();
		super.dispose();
		}
		
	private static native long _open(final String s,final String m);
	private static native void _close(final long ptr);
	private static native int _getline(final long ptr,int delim,long ks);

}
