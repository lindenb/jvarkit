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
package com.github.lindenb.jvarkit.setfile;

import java.io.PrintWriter;
import java.util.AbstractList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.util.Locatable;

/* 'SetFile' for rvtest */
public interface SetFileRecord extends List<Locatable> {
	public static final String FILE_EXTENSION=".setfile";
	public String getName();
	
	/** return all chromosomes in this setFileRecord */
	public default Set<String> getChromosomes() {
		return this.stream().map(X->X.getContig()).collect(Collectors.toSet());
	}
	
	/** return all sum of length on reference */
	public default long getLongSumOfLengthOnReference() {
		return this.stream().mapToLong(R->R.getLengthOnReference()).sum();
	}
	/** return true if any item intersect with 'loc' */
	public default boolean overlaps(final Locatable loc) {
		return this.stream().anyMatch(R->R.overlaps(loc));
	}
	
	public default void println(final PrintWriter pw) {
		pw.write(getName());
		for(int i=0;i< size();i++) {
			final Locatable rec = get(i);
			pw.write(i==0?"\t":",");
			pw.write(rec.getContig());
			pw.write(":");
			pw.write(String.valueOf(rec.getStart()));
			pw.write("-");
			pw.write(String.valueOf(rec.getEnd()));
			}
		pw.write("\n");
	}
	
	public static SetFileRecord create(final String name,final List<Locatable> intervals) {
		return new SetRecordImpl(name,intervals);
		}
	
	
	class SetRecordImpl extends AbstractList<Locatable>
		implements SetFileRecord {
		private final String name;
		private final List<Locatable> intervals;
		SetRecordImpl(final String name,final List<Locatable> intervals) {
			this.name = name;
			this.intervals = intervals;
			}
		@Override
		public String getName() {
			return this.name;
			}
		@Override
		public int size() {
			return this.intervals.size();
			}
		@Override
		public Locatable get(int i) {
			return this.intervals.get(i);
			}
		
		@Override
		public boolean equals(Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof SetFileRecord)) return false;
			final SetFileRecord rec = SetFileRecord.class.cast(obj);
			if(this.size()!=rec.size()) return false;
			if(!this.getName().equals(rec.getName())) return false;
			for(int i=0;i< this.size();i++) {
				final Locatable a = this.get(i);
				final Locatable b = rec.get(i);
				if(!a.contigsMatch(b)) return false;
				if(a.getStart()!=b.getStart()) return false;
				if(a.getEnd()!=b.getEnd()) return false;
				}
			return true;
			}
		
		@Override
		public int hashCode() {
			return this.name.hashCode() + this.intervals.hashCode();
			}
		
		@Override
		public String toString() {
			return getName()+"\t"+ this.stream().
					map(R->R.getContig()+":"+R.getStart()+"-"+R.getEnd()).
					collect(Collectors.joining(","));
			}
		}
	}
