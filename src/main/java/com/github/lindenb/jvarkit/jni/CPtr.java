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
package com.github.lindenb.jvarkit.jni;

public class CPtr {
	private long ptr = 0;
	
	public CPtr()
		{
		this(0L);
		}
	
	public CPtr(long n) {
		this.ptr = n;
		}
	
	@Override
	protected void finalize() throws Throwable {
		dispose();
		super.finalize();
		}
	
	public void dispose() {
		this.set(0L);
		}
	
	public void assertNotNull() {
		if(isNull()) throw new NullPointerException("JNI Pointer is null !");
	}
	
	public long getMustBeNotNull() {
		assertNotNull();
		return this.get();
		}

	
	public long get() {
		return this.ptr;
		}
	
	public void set(long n) {
		this.ptr = n;
	}
	
	public void setNull() {
		set(0L);
	}
	
	public boolean isNull() {
		return this.ptr<=0L;
		}
	@Override
	public int hashCode() {
		return Long.hashCode(this.ptr);
		}
	
	@Override
	public String toString() {
		return "ptr "+(isNull()?"(null)":" adr:"+this.get());
		}
 	}
