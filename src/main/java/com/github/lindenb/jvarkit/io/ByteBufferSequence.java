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
package com.github.lindenb.jvarkit.io;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.UnsupportedEncodingException;
import java.util.Arrays;
import java.util.function.IntUnaryOperator;

import com.github.lindenb.jvarkit.lang.SubSequence;

import htsjdk.samtools.Defaults;
/**
 * A mustable string, using bytes instead of char
 * @author lindenb
 *
 */
public class ByteBufferSequence extends OutputStream implements Cloneable,CharSequence,Appendable,Comparable<ByteBufferSequence> {
	protected byte[] buf;
	protected int count;
	public ByteBufferSequence() {
		this(Defaults.NON_ZERO_BUFFER_SIZE);
		}
	public ByteBufferSequence(int capacity) {
		buf=new byte[capacity];
		count =0;
		}
	
	public ByteBufferSequence clear() {
		count=0;
		return this;
		}
	
	@Override
	protected ByteBufferSequence clone() {
		return new ByteBufferSequence(this.buf,0,this.size());
		}
	
	
	public ByteBufferSequence(final CharSequence seq) {
		this(seq.length());
		if(seq instanceof ByteBufferSequence) {
			final ByteBufferSequence o = ByteBufferSequence.class.cast(seq);
			System.arraycopy(o.buf,0, this.buf,0, o.count);
			}
		else
			{
			for(int i=0;i< seq.length();i++) {
				write(seq.charAt(i));
				}
			}
		}
	
	public ByteBufferSequence(byte[] array,int offset,int length) {
		this(length);
		System.arraycopy(array, offset, this.buf, 0, length);
		this.count=length;
		}
	
	public ByteBufferSequence(byte[] array) {
		this(array,0,array.length);
		}
	
	 public void ensureCapacity(int minCapacity) {
	    if (minCapacity > buf.length) {
	        buf = Arrays.copyOf(buf,Math.max(minCapacity, (int)( buf.length*1.5)));
	    	}
	    }

	 public void transform(IntUnaryOperator fun) {
		 for(int i=0;i< this.count;++i) {
			 buf[i] = (byte)fun.applyAsInt(buf[i]);
		 	}
	 	}
	 
	 @Override
	 public void write(int b) {
        ensureCapacity(count + 1);
        buf[count] = (byte) b;
        count += 1;
	    }
	
	@Override
	public void write(byte[] b, int off, int len) throws IOException {
		ensureCapacity(count + len);
		System.arraycopy(b, off, this.buf, count, len);
		count+=len;
		}
	 
	public int size() {
		return count;
		}
	
	
	@Override
	public final int length() {
		return this.size();
	}

	@Override
	public char charAt(int index) {
		return (char)this.buf[index];
	}
	
	public void setByteAt(int i,byte b) {
		if(i<0 ||i>=count) throw new IllegalArgumentException();
		buf[i]=b;
		}
	
	public void setCharAt(int i,char c) {
		setByteAt(i,(byte)c);
		}

	@Override
	public CharSequence subSequence(int start, int end) {
		return new SubSequence(this,start,end);
	}

	@Override
	public ByteBufferSequence append(final CharSequence csq) throws IOException {
		return append(csq,0,csq.length());
		}

	@Override
	public ByteBufferSequence append(CharSequence csq, int start, int end) throws IOException {
		for(int i=start;i< end;i++) {
			append(csq.charAt(i));
			}
		return this;
		}

	@Override
	public ByteBufferSequence append(char c) throws IOException {
		write((char)c);
		return this;
		}
	
	@Override
	public int hashCode() {
		int result = 1;
        for(int i=0;i< size();i++) {
        	byte element = this.buf[i];
            result = 31 * result + element;
        	}
        return result;
		}
	
	@Override
	public boolean equals(Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof ByteBufferSequence)) return false;
		final ByteBufferSequence o=ByteBufferSequence.class.cast(obj);
	    if(this.size()!=o.size()) return false;
	    for(int i=0;i< this.count;i++) {
            byte c1 = this.buf[i];
            byte c2 = o.buf[i];
            if (c1 != c2) return false;
	  		}
        return true;
		}
	
	@Override
	public int compareTo(final ByteBufferSequence o) {
	  final int lim = Math.min(this.size(), o.size());
	  for(int i=0;i< lim;i++) {
            byte c1 = this.buf[i];
            byte c2 = o.buf[i];
            if (c1 != c2) {
                return c1 - c2;
            	}
	  		}
       return size() - o.size();	
       }
	
	 public byte[] toByteArray() {
        return Arrays.copyOf(buf, count);
    	}
	 
	 public ByteArrayInputStream toByteArrayInputStream() {
		 return new ByteArrayInputStream(this.buf, 0, this.size());
	 	}
	 
	 public boolean startsWith(String s) {
		 if(s.length()> this.length()) return false;
		 for(int i=0;i< s.length();i++) {
			if(s.charAt(i)!=this.charAt(i))  return false; 
		 	}
		 return true;
	 	}
	
	@Override
	public String toString() {
		try {
			return new String(toByteArray(),"UTF-8");
		} catch (UnsupportedEncodingException err) {
			throw new RuntimeException(err);
			}
		}
	}
