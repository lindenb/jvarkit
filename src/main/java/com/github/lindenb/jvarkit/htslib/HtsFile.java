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
package com.github.lindenb.jvarkit.htslib;

import java.io.Closeable;
import java.io.IOException;
import java.nio.file.Path;
import htslib.Htslib;
import com.github.lindenb.jvarkit.jni.CPtr;


public abstract class HtsFile extends CPtr implements Closeable {
	private final String filename;
	public HtsFile(final Path s,final String m) throws IOException {
		this(s.toString(),m);
		}
	
	public HtsFile(final String s,final String m) throws IOException {
		super(Htslib.hts_hopen(s,m));
		this.filename = s;
		if(isNull()) throw new IOException("Cannot open "+s);
		}
	
	public boolean isOpen() {
		return !super.isNull();
		}
	
	@Override
	public final void close() {
		this.dispose();
		}
		
	@Override
	public void dispose() {
		if(isOpen()) Htslib.hts_hclose(this.getPtr());
		setNull();
		super.dispose();
		}
	
	public String getFilename() {
		return filename;
		}
	
	@Override
	public String toString() {
		return this.getClass().getName()+"(" + getFilename() + ")";
		}
	}
