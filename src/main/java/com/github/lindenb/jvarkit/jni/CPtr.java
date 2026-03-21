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
package com.github.lindenb.jvarkit.jni;

/** wrapper for C pointers */
public interface CPtr {
	
	/** set ptr to 0L, should be the place to free memory if needed */
	public void disposePtr();
	
	
	public default long getPtrMustBeNotNull() {
		long n = getPtr();
		if(n==0L) throw new IllegalArgumentException("Ptr is null");
		return this.getPtr();
		}
	
	public long getPtr();
	public void setPtr(final long n);
	
	public default void setPtrToNull() {
		setPtr(0L);
		}
	
	public default boolean isPtrNull() {
		return getPtr()==0L;
		}
	public default boolean isPtrNotNull() {
		return !isPtrNull();
		}
 	}
