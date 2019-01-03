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
package com.github.lindenb.jvarkit.util.vcf;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.ParameterException;

/**
 * 
 * Class wrapping chrom/pos
 *
 */
public class ContigPos implements htsjdk.samtools.util.Locatable,Comparable<ContigPos>{
	
	/** converter for JCommander */
	public static class Converter implements IStringConverter<ContigPos>
		{
		@Override
		public ContigPos convert(final String s)
			{
			if(s==null) return null;
			try {
				return new ContigPos(s.trim());
				}
			catch(Exception err) {
				throw new ParameterException("Cannot convert "+s+" to Contig/pos",err);
				}
			}
		}
	

	private final String contig;
	private final int pos;
	public ContigPos(final String contig,final int pos)
		{
		this.contig = contig;
		this.pos = pos;
		}

	/** parse using the syntax : 'chr:pos' */
	public ContigPos(final String str)
		{
		int colon = str.indexOf(":");
		if( colon == -1) throw new IllegalArgumentException("no colon in "+str);
		if( colon == 0) throw new IllegalArgumentException("no contig in \""+str+"\"");
		this.contig = str.substring(0,colon);
		this.pos = Integer.parseInt(str.substring(colon+1));
		}
	
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + contig.hashCode();
		result = prime * result + pos;
		return result;
	}

	@Override
	public int compareTo(final ContigPos o) {
		int i = this.contig.compareTo(o.contig);
		if( i!=0 ) return i; 
		return this.pos - o.pos;
		}
	
	@Override
	public boolean equals(final Object obj) {
		if (this == obj) {
			return true;
		}
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		final ContigPos other = (ContigPos) obj;
		return  this.pos==other.pos &&
				this.contig.equals(other.contig)
				;
	}

	@Override
	public String getContig() {
		return contig;
		}
	public int getPos() {
		return pos;
		}
	/** same as getPos  */
	@Override
	public int getStart() {
		return getPos();
		}
	/** getEnd()==getStart()  */
	@Override
	public int getEnd() {
		return getPos();
		}
	@Override
	public String toString() {
		return getContig()+":"+getPos();
		}
	}
