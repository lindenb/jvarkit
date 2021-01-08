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
package com.github.lindenb.jvarkit.jni;

/** wrapper for C pointers */
public class CPtr {
	private long ptr = 0L;
	
	public CPtr()
		{
		this(0L);
		}
	
	public CPtr(final long n) {
		this.ptr = n;
		}
	
	@Override
	protected void finalize() throws Throwable {
		dispose();
		super.finalize();
		}
	
	/** set ptr to 0L, should be the place to free memory if needed */
	public void dispose() {
		this.setPtr(0L);
		}
	
	public void assertNotNull() {
		if(isNull()) throw new NullPointerException("JNI Pointer is null !");
	}
	
	public long getMustBeNotNull() {
		assertNotNull();
		return this.getPtr();
		}

	
	public long getPtr() {
		return this.ptr;
		}
	
	public void setPtr(final long n) {
		this.ptr = n;
	}
	
	public void setNull() {
		setPtr(0L);
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
		return String.valueOf(this.getClass().getName())+"("+(isNull()?"null":String.format("0x%08X",this.getPtr())+")");
		}
 	}
