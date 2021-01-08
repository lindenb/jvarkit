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
package com.github.lindenb.jvarkit.samtools;

import java.io.IOException;
import java.io.StringWriter;
import java.io.Writer;
import java.util.function.Function;

import com.beust.jcommander.IStringConverter;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqConstants;
import htsjdk.samtools.util.RuntimeIOException;

/**
 * How to print a Read using a printf-like function
 */
public interface SAMRecordNaming extends Function<SAMRecord, String> {
	
	public static final String OPT_DESC="How to print a Read. A format 'a la C-printf'. "
			+ "%% :% , %n:read name, %s: read bases, %q: read quals, %f : read flags,"
			+ "%m: mapq, %c: contig, %b: start, %B: unclipped start, %e: end, %E: unclipped end,"
			+ "%I: read group id, %N: sample name,%S: SAM String.";

	/** jcommander converter */
	public static class StringConverter implements IStringConverter<SAMRecordNaming> {
		@Override
		public SAMRecordNaming convert(final String s) {
			return SAMRecordNaming.compile(s);
			}
		}
	
	/** compile a new  SAMRecordNaming */
	public static SAMRecordNaming compile(final String s) {
		return new SAMRecordNamingImpl(s);
	}
	
	/** print samrecord to w using the current pattern */ 
	public void print(final Writer w,final SAMRecord t);
	
	@Override
	public default String apply(final SAMRecord t) {
		final StringWriter sw = new StringWriter();
		print(sw,t);
		return sw.toString();
		}
	
	static class SAMRecordNamingImpl implements SAMRecordNaming {
		final String pattern;
		SAMRecordNamingImpl(final String s) {
			this.pattern = (s==null?"":s);
		}
		
	
		
	@Override
	public void print(final Writer w,final SAMRecord rec) {
		try  {
			if(rec==null || this.pattern==null) return;
			int i=0;
			while(i< this.pattern.length()) {
				if(this.pattern.charAt(i)=='%') {
					if(i+1==this.pattern.length()) return;
					i++;
					switch(this.pattern.charAt(i))
						{
						case '%': w.append('%');break;
						case 'n': {
							w.append(rec.getReadName());
							if(rec.getReadPairedFlag()) {
					        	if(rec.getFirstOfPairFlag()) w.append(FastqConstants.FIRST_OF_PAIR);
					        	if(rec.getSecondOfPairFlag()) w.append(FastqConstants.SECOND_OF_PAIR);
							}
							break;
						}
						case 's': w.append(rec.getReadString());break;
						case 'q': w.append(rec.getBaseQualityString());break;
						case 'f': w.append(String.valueOf(rec.getFlags()));break;
						case 'm': w.append(String.valueOf(rec.getMappingQuality()));break;
						case 'c': if(!rec.getReadUnmappedFlag()) w.append(rec.getContig());break;
						case 'b': if(!rec.getReadUnmappedFlag()) w.append(String.valueOf(rec.getStart()));break;
						case 'B': if(!rec.getReadUnmappedFlag()) w.append(String.valueOf(rec.getUnclippedStart()));break;
						case 'e': if(!rec.getReadUnmappedFlag()) w.append(String.valueOf(rec.getEnd()));break;
						case 'E': if(!rec.getReadUnmappedFlag()) w.append(String.valueOf(rec.getUnclippedEnd()));break;
						case 'I': {SAMReadGroupRecord rg = rec.getReadGroup();if(rg!=null)  w.append(rg.getId());break;}
						case 'N': {SAMReadGroupRecord rg = rec.getReadGroup();if(rg!=null && rg.getSample()!=null)  w.append(rg.getSample());break;}
						case 'A':  w.append(rec.getSAMString());break;
						default:break;
						}
					i++;
					}
				else if(this.pattern.charAt(i)=='\\') {
					if(i+1==this.pattern.length()) return;
					i++;
					switch(this.pattern.charAt(i))
						{
						case '\\': w.append('\\');break;
						case '\'': w.append('\'');break;
						case '\"': w.append('\"');break;
						case 't': w.append('\t');break;
						case 'n': w.append('\n');break;
						default:break;
						}
					i++;
					}
				else
					{
					w.append(this.pattern.charAt(i));
					i++;
					}
				}
			} catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	
	@Override
	public boolean equals(final Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof SAMRecordNamingImpl)) return false;
		return this.pattern.equals(SAMRecordNamingImpl.class.cast(obj).pattern);
		}
	
	@Override
	public int hashCode() {
		return this.pattern.hashCode();
		}
	
	@Override
	public String toString() {
		return this.pattern;
		}
	}
}
