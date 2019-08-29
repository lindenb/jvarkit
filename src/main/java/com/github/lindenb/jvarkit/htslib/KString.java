package com.github.lindenb.jvarkit.htslib;

import com.github.lindenb.jvarkit.jni.CPtr;

public class KString
extends CPtr  
	implements CharSequence,
	Comparable<CharSequence>
{
public KString() {
	super(Htslib.bind_kstring_create());
	}
	
@Override
public char charAt(int index) { return (char)Htslib.bind_kstring_at(this.get(),index);}
@Override
public int length() { return Htslib.bind_kstring_len(this.get());}
@Override
public CharSequence subSequence(int start, int end) {
	int L=end-start;
	final StringBuilder sb=new StringBuilder(L);
	for(int i=0;i< L;i++) sb.append(charAt(start+i));
	return sb;
	}
@Override
public String toString() {
	return subSequence(0,length()).toString();
	}

@Override
public void dispose() {
	if(!isNull()) Htslib.bind_kstring_release(this.get());
	super.dispose();
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

}
