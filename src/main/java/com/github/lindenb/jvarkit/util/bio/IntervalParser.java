/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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


History:
* 2017 creation

*/
package com.github.lindenb.jvarkit.util.bio;

import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.Set;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;

public class IntervalParser {
	private SAMSequenceDictionary dict=null;
	private boolean raiseExceptionOnError=true;
	private boolean tryToFixContigName=true;
	
	public IntervalParser()
		{
		this(null);
		}
	public IntervalParser(final SAMSequenceDictionary dict)
		{
		this.dict = dict;
		}

	/** there is a dictionary and we cannot find the SAMSequenceRecord
	 * try to find an equivalent chromosome e.g: 1->chr1 , chr2 -> 2, MT->chrM
	 * @param tryToFixContigName
	 * @return
	 */
	public IntervalParser setFixContigName(boolean tryToFixContigName) {
		this.tryToFixContigName = tryToFixContigName;
		return this;
		}
	
	public boolean isFixContigName() {
		return tryToFixContigName;
		}
	
	public IntervalParser setDictionary(final SAMSequenceDictionary dict) {
		this.dict = dict;
		return this;
		}
	public boolean hasDictionary() {
		return getDictionary()!=null;
	}
	public SAMSequenceDictionary getDictionary() {
		return this.dict;
		}
	public IntervalParser setRaiseExceptionOnError(boolean raiseExceptionOnError) {
		this.raiseExceptionOnError = raiseExceptionOnError;
		return this;
		}
	
	private Interval returnErrorOrNullInterval(final String message) {
		if(raiseExceptionOnError)
			{
			throw new IllegalArgumentException(message);
			}
		return null;
		}
	
	private Set<String> findAlternativeContigNames(String contig)
		{
		if(!this.tryToFixContigName) return Collections.emptySet();
		final Set<String> set = new LinkedHashSet<>();
		if(contig.startsWith("chr"))
			{
			set.add(contig.substring(3));
			}
		else if(!contig.startsWith("chr"))
			{
			set.add("chr"+contig);
			}
			
		if(contig.equals("chrM")) set.add("MT");
		if(contig.equals("MT")) set.add("chrM");
		
		return set;
		}
	
	private SAMSequenceRecord getSAMSequenceRecord(final String contig)
		{
		if(!hasDictionary()) return null;
		SAMSequenceRecord rec= getDictionary().getSequence(contig);
		if(rec!=null) return rec;
		//case problem ?
		if(this.isFixContigName() )
			{
			for(final SAMSequenceRecord ssr: getDictionary().getSequences())
				{
				if(ssr.getSequenceName().equalsIgnoreCase(contig)) return ssr;
				}
			}
		for(final String c: findAlternativeContigNames(contig)) {
			rec= getDictionary().getSequence(c);
			if(rec!=null) return rec;
			}
		return null;
		}
	
	public Interval parse(final String s)
		{
		final int colon=s.indexOf(':');
		if(colon<1 || colon+1==s.length()) {
			return returnErrorOrNullInterval("Cannot find colon in "+s);
			}
		final int hyphen=s.indexOf('-',colon+1);
		if(hyphen==-1) return returnErrorOrNullInterval("Cannot find hyphen in "+s);
			
		int start,end;
		try {
			SAMSequenceRecord ssr=null;
			String chrom=s.substring(0,colon);
			if(hasDictionary())
				{
				ssr= getSAMSequenceRecord(chrom);
				if(ssr==null) return returnErrorOrNullInterval("Cannot find chromosome "+chrom+" in dictionary.");
				chrom = ssr.getSequenceName();
				}
			start=Integer.parseInt(s.substring(colon+1,hyphen));
			end=Integer.parseInt(s.substring(hyphen+1));
			if(end<start) return null;
			return new Interval(s.substring(0,colon), start, end,false,s);
		} catch (final Exception e)
			{
			return null;
			}
		}
	
	}
