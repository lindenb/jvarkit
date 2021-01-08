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
package com.github.lindenb.jvarkit.lang;


import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.OptionalDouble;
import java.util.OptionalInt;
import java.util.OptionalLong;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import com.github.lindenb.jvarkit.io.IOUtils;

import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;

public class StringUtils extends StringUtil {
	private static final DecimalFormat NICE_INT_FORMAT = new DecimalFormat("###,###");


private static <T> boolean _is(final String s,Function<String,T> fun,final Predicate<T> validator) {
	if(isBlank(s)) return false;
	try { T t=fun.apply(s);return t!=null && (validator==null || validator.test(t));}
	catch(final Throwable err) {
		return false;
		}
	}
	
public static boolean isDouble(final String s) {
	return parseDouble(s).isPresent();
	}
public static boolean isFloat(final String s) {
	return _is(s,Float::parseFloat,null);
	}
public static boolean isInteger(final String s) {
	return parseInt(s).isPresent();
	}
public static boolean isLong(final String s) {
	return parseLong(s).isPresent();
	}

public static OptionalInt parseInt(final String s) {
	if(StringUtils.isBlank(s)) return OptionalInt.empty();
	try {
		return OptionalInt.of(Integer.parseInt(s));
	} catch(final NumberFormatException err) {
		return OptionalInt.empty();
	}
}

public static OptionalLong parseLong(final String s) {
	if(StringUtils.isBlank(s)) return OptionalLong.empty();
	try {
		return OptionalLong.of(Long.parseLong(s));
	} catch(final NumberFormatException err) {
		return OptionalLong.empty();
	}
}

public static OptionalDouble parseDouble(final String s) {
	if(StringUtils.isBlank(s)) return OptionalDouble.empty();
	try {
		return OptionalDouble.of(Double.parseDouble(s));
	} catch(final NumberFormatException err) {
		return OptionalDouble.empty();
	}
}


/**
 * returns a string that is the rest of a given string after a given substring.
 * @param haystack The string to be evaluated. Part of this string will be returned.
 * @param needle The substring to search for. Everything after the first occurence ofneedle inhaystack will be returned.
 * @return a string, empty string if the needle was not found in haystack
 */
public static String substringAfter(final String haystack ,final String needle ) {
	final int i = haystack.indexOf(needle);
	if(i==-1) return "";
	return haystack.substring(i+needle.length());
	}
/**
 * returns a string that is the rest of a given string before a given substring.
 * @param haystack The string to be evaluated. Part of this string will be returned.
 * @param needle The substring to search for. Everything before the first occurance ofneedle inhaystack will be returned.
 * @return a string, empty string if the needle was not found in haystack
 */
public static String substringBefore(final String haystack ,final String needle ) {
	final int i = haystack.indexOf(needle);
	if(i==-1) return "";
	return haystack.substring(0,i);
	}
/** append line numbers to code */
public static String beautifyCode(final String sourceCode)
	{
	final String codeLines[] = sourceCode.split("[\n]");
	final StringBuilder codeWithLineNumber = new StringBuilder(sourceCode.length()+(15*codeLines.length));
	for(int nLine=0;nLine < codeLines.length;++nLine)
		{
		codeWithLineNumber.
			append(nLine==0?"":"\n").
			append(String.format("%10d  ",(nLine+1))+codeLines[nLine])
			;
		}
	return codeWithLineNumber.toString();
	}
/** repeat  char `c` for  `n` time */
public static String repeat(final int n,final char c) {
	if(n<0) throw new IllegalArgumentException("repeat: n is negative");
	return repeatCharNTimes(c, n);
	}
/** escape Http using UTF-8 */
public static String escapeHttp(final String str) {
	if(str==null) return null;
	try {
		return URLEncoder.encode(str,"UTF-8");
		}
	catch(final UnsupportedEncodingException err) {
		throw new IllegalArgumentException(err);
		}
	}
/** return first letter as UpperCase and reminder to lowerCase . Return null if argument is null */
public static String toTitle(final String s) {
	if(s==null) return null;
	if(s.isEmpty()) return s;
	if(s.length()==1) return s.toUpperCase();
	return s.substring(0, 1).toUpperCase() + s.substring(1).toLowerCase();
	}
/** escape C string , returns null is str is null */
public static String escapeC(final String str) {
	if(str==null) return null;
	
	int i=0;
	while(i< str.length())
		{
		char c = str.charAt(i);
		switch(c)
			{
			case '\'': 
			case '\"':
			case '\n':
			case '\t':
			case '\\':
				{
				final StringBuilder sb = new StringBuilder(str.length()+1);
				sb.append(str.substring(0, i));
				while(i< str.length())
					{
					c = str.charAt(i);
					switch(c)
						{
						case '\'': sb.append("\\\'");break;
						case '\"': sb.append("\\\"");break;
						case '\n': sb.append("\\n");break;
						case '\t': sb.append("\\t");break;
						case '\\': sb.append("\\\\");break;
						default: sb.append(c);break;
						}
					i++;
					}
				return sb.toString();
				}
			default: break;
			}
		i++;
		}
	return str;
	}

/** return long number with comma as thousand separator */
public static final String niceInt(long n) {
	return StringUtils.NICE_INT_FORMAT.format(n);
	}

/** return left part of a string */
public static final String left(final CharSequence s,int len) {
	if(s.length()<=len) return s.toString();
	return s.subSequence(0, len).toString();
	}
/** return right part of a string */
public static final String right(final CharSequence s,int len) {
	if(s.length()<=len) return s.toString();
	return s.subSequence(s.length()-len,s.length()).toString();
	}
/** return md5 of a string */
public static final String md5(final String in) {
	final MessageDigest _md5;
	try {
		_md5 = java.security.MessageDigest.getInstance("MD5");
	} catch (final NoSuchAlgorithmException e) {
		throw new RuntimeException("MD5 algorithm not found", e);
		}
	
	_md5.reset();
	_md5.update(in.getBytes());
	String s = new java.math.BigInteger(1, _md5.digest()).toString(16);
	if (s.length() != 32) {
		final String zeros = "00000000000000000000000000000000";
		s = zeros.substring(0, 32 - s.length()) + s;
		}
	return s;
	}

/** return true if s ends with any of the suffixes */
public static final boolean endsWith(final String s,String...suffixes) {
	if(s==null) throw new IllegalArgumentException("s is null");
	for(final String suff:suffixes) if(s.endsWith(suff)) return true;
	return false;
	}

/** unescape quoted C string. Starting and trailing quote have already been removed */
public static String unescapeC(final String s)
	{
	if(s==null) throw new IllegalArgumentException("s is null");
	final StringBuilder b=new StringBuilder(s.length());
	int i=0;
	while(i<s.length())
		{
		if(s.charAt(i)=='\\')
			{
			if( i+1== s.length())  break;
			++i;
			switch(s.charAt(i))
				{
				case 'n': b.append("\n");break;
				case 'r': b.append("\r");break;
				case 't': b.append("\t");break;
				case '\\': b.append("\\");break;
				case '\'': b.append("\'");break;
				case '\"': b.append("\"");break;
				default: break;//ignore
				}
			}
		else
			{
			b.append(s.charAt(i));
			}
		++i;
		}
	return b.toString();
}

/** convert duration to a nice string */
public static String niceDuration(long durationMillisec) {
	long n =durationMillisec/1000;
	if(n<60)
		{
		return n+" second"+(n<2?"":"s");
		}
	n/=60;//minutes
	if(n< 60)
		{
		return n+" minute"+(n<2?"":"s");
		}
	n/=60;//hours
	if(n< 24)
		{
		return n+" hour"+(n<2?"":"s");
		}
	n/=24;
	
	if(n< 365)
		{
		return n+" day"+(n<2?"":"s");
		}
	n/=365;
	
	return n+" year"+(n<2?"":"s");
	}
/** convert file size to a nice string 
 * source: https://stackoverflow.com/questions/3758606 */
public static String niceFileSize(final long bytes) {
	final boolean si = false;
	final  int unit = si ? 1000 : 1024;
    if (bytes < unit) return bytes + " B";
    final int exp = (int) (Math.log(bytes) / Math.log(unit));
    final String suffixes = (si ? "kMGTPE" : "KMGTPE");
    if(exp-1 >= suffixes.length()) return String.valueOf(bytes);
    final String pre = suffixes.charAt(exp-1) + (si ? "" : "i");
    return String.format("%.1f %sB", bytes / Math.pow(unit, exp), pre);
	}

public static String ifBlank(final String s,final String def) {
	return isBlank(s)?def:s;
	}

/** if s is null return empty string. replace all CR per space. Trim and Replace all consecutive spaces */
public static String normalizeSpaces(final String s) {
	if(s==null) return "";
	return s.replaceAll("[\\s]+", " ").trim();
	}
/** return date using format yyyyMMdd_HHmmss */
public static String now() {
	return new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
	}

/** return index of 'needle' in haystack, case insensitive */
public static int indexOfIgnoreCase(final String haystack,final String needle) {
	//return default behavior if any empty 
	if (needle.isEmpty() || haystack.isEmpty()) return haystack.indexOf(needle);
	
	for (int i = 0; i + needle.length() <= haystack.length(); ++i) {
		int j=0;
		for(j=0;
				j< needle.length() && 
				Character.toLowerCase(haystack.charAt(i+j)) == 
				Character.toLowerCase(needle.charAt(j));
			j++) {
			//nothing
			}
		if(j==needle.length()) return i;
	 	}
    return -1;
	}

/** compress string as gzipped bytes */
public byte[] compressString(final String str) {
	try(ByteArrayOutputStream baos = new ByteArrayOutputStream(Math.max(32,str.length()/2))) {
		try(GZIPOutputStream gzout = new GZIPOutputStream(baos)) {
			gzout.write(str.getBytes());
			gzout.finish();
			gzout.flush();
			}
		return baos.toByteArray();
		}
	catch(final IOException err) {
		throw new RuntimeIOException(err);
		}
	}
/** compress string as gzipped bytes */
public String uncompressString(final byte[] compressed) {
	try(ByteArrayInputStream bis = new ByteArrayInputStream(compressed)) {
		try(GZIPInputStream gzin = new GZIPInputStream(bis)) {
			try(ByteArrayOutputStream bos = new ByteArrayOutputStream()) {
				IOUtils.copyTo(gzin, bos);
				return new String(bos.toByteArray());
				}
			}
		}
	catch(final IOException err) {
		throw new RuntimeIOException(err);
		}
	}

}