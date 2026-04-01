package com.github.lindenb.jvarkit.htslib;

import java.io.IOException;

import com.github.lindenb.jvarkit.jni.AbstractCPtr;

public class KString extends AbstractCPtr implements CharSequence, Comparable<CharSequence>, Appendable
	{
	public KString() {
		this(KString.kstring_create(),true);
		}
	
	public KString(long ptr,boolean managed) {
		super(ptr,managed);
		}
		
	@Override
	public char charAt(int index) { return (char)kstring_at(this.getPtrMustBeNotNull(),index);}
	@Override
	public int length() { return KString.kstring_len(this.getPtrMustBeNotNull());}
	
	@Override
	public CharSequence subSequence(int start, int end) {
		final int L=end-start;
		final StringBuilder sb=new StringBuilder(L);
		for(int i=0;i< L;i++) sb.append(charAt(start+i));
		return sb;
		}
	
	@Override
	public String toString() {
		return subSequence(0,length()).toString();
		}
	
	@Override
	public void disposePtr() {
		if(!isPtrNull()) KString.kstring_release(this.getPtr());
		setPtrToNull();
		}
	
	public int compareTo(final CharSequence anotherString) {
		final int len1 = this.length();
		final int len2 = anotherString.length();
		final int lim = Math.min(len1, len2);
	
	    int k = 0;
	    while (k < lim) {
	        final char c1 = this.charAt(k);
	        final char c2 = anotherString.charAt(k);
	        if (c1 != c2) {
	            return c1 - c2;
	        }
	        k++;
	    }
	    return len1 - len2;
		}
	
	public boolean isEmpty() {
		return length()>0;
		}
	public int hashCode() {
	    int h = 0;
	    final int L=this.length();
	     for (int i = 0; i < L; i++) {
	            h = 31 * h + (int)charAt(i);
	        }
	    return h;
		}
	
	@Override
	public Appendable append(char c) throws IOException {
		kstring_putc(getPtrMustBeNotNull(),(byte)c);
		return this;
		}
	@Override
	public Appendable append(CharSequence csq) throws IOException {
		return append(csq,0,csq.length());
		}
	@Override
	public Appendable append(CharSequence csq, int start, int end) throws IOException {
		while(start < end) {
			append(csq.charAt(start++));
			}
		return this;
		}
	
	
	private static native long kstring_create();
	private static native void kstring_release(long ptr);
	private static native int kstring_len(long ptr);
	private static native byte kstring_at(long ptr,int idx);
	private static native byte kstring_putc(long ptr,byte c);
	}
