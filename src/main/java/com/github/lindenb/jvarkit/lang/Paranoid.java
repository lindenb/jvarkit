/*

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
The MIT License (MIT)
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
package com.github.lindenb.jvarkit.lang;

import java.util.Arrays;
import java.util.stream.Collectors;

/** on-time assertion checker */
public interface Paranoid {

public default Paranoid assertTrue(boolean b) {
	return assertTrue(b,"True Assertion failed");
	}
public default Paranoid assertTrue(boolean b,final String msg) {
	if(!b) raise(msg);
	return this;
	}

public default Paranoid assertFalse(boolean b) {
	return assertFalse(b,"FALSE Assertion failed");
	}
public default Paranoid assertFalse(boolean b,final String msg) {
	if(b) raise(msg);
	return this;
	}
public default Paranoid assertLt(final int a,final int b) {
	if(!(a<b)) raise("assertion "+ a+" < "+ b +" failed");
	return this;
	}
public default Paranoid assertLe(final int a,final int b) {
	if(!(a<=b)) raise("assertion "+ a+" <= "+ b +" failed");
	return this;
	}
public default Paranoid assertGt(final int a,final int b) {
	if(!(a>b)) raise("assertion "+ a+" > "+ b +" failed");
	return this;
	}
public default Paranoid assertGe(final int a,final int b) {
	if(!(a>=b)) raise("assertion "+ a+" >= "+ b +" failed");
	return this;
	}

public default Paranoid assertEq(final int a,final int b) {
	if(!(a==b)) raise("assertion "+ a+" == "+ b +" failed");
	return this;
	}
public default Paranoid assertNe(final int a,final int b) {
	if(a==b) raise("assertion "+ a+" != "+ b +" failed");
	return this;
	}


public Paranoid raise(final String msg);

public default void dump() {
}


static abstract class AbstractParanoid implements Paranoid {
	
	
	public String where() {
		String msg="undefined";
		try {
			final StackTraceElement array[]=Thread.currentThread().getStackTrace();
			if(array!=null && array.length>0)
				{
				final String s= Arrays.stream(array).
					filter(STE->!STE.getClassName().startsWith(Paranoid.class.getName())).
					//filter(STE->STE.getClassName().startsWith("com.github.lindenb.")).
					map(STE->STE.getClassName()+" in \""+STE.getFileName()+"\" "+STE.getMethodName()+" line "+STE.getLineNumber()).
					collect(Collectors.joining("\n"));
					;
				return s.isEmpty()?msg:s;
				}
			}
		catch(final Throwable err) {
			}
		return msg;
		}
	}

static class Silent extends AbstractParanoid {
	@Override
	public Paranoid raise(String msg) {
		return this;
		}
	}
static class Eager extends AbstractParanoid {
	@Override
	public Paranoid raise(String msg) {
		throw new AssertionError(String.valueOf(msg)+" "+where());
		}
	}

public static Paranoid createSilentInstance() { 
	return new Silent();
	}

public static Paranoid createThrowingInstance() { 
	return new Eager();
	}
}
