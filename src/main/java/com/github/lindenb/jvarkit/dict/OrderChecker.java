/*
The MIT License (MIT)

Copyright (c) 2022 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.dict;

import java.util.HashSet;
import java.util.Set;
import java.util.function.UnaryOperator;

import com.github.lindenb.jvarkit.lang.JvarkitException;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;

public class OrderChecker<T extends Locatable> implements UnaryOperator<T>{
private final Set<String> contig_seen;
private final SAMSequenceDictionary dict;
private final boolean strictDictOrder;
private String prev_contig = null;
private int prev_start = -1;

public OrderChecker() {
	this(null,false);
}

public OrderChecker(final SAMSequenceDictionary dict) {
	this(dict,true);
}

public OrderChecker(final SAMSequenceDictionary dict,boolean strictDictOrder) {
	this.dict = (dict.isEmpty()?null:dict);
	this.strictDictOrder = strictDictOrder;
	if(this.strictDictOrder) {
		if (this.dict==null) throw new IllegalArgumentException("strict ordering but dict is null or empty.");
		this.contig_seen = null;
		}
	else {
		this.contig_seen = new HashSet<>();
		}
	}

private int getTid(final String c) {
	final SAMSequenceRecord ssr = dict.getSequence(c);
	if(ssr==null) throw new JvarkitException.ContigNotFoundInDictionary(c, this.dict);
	return ssr.getSequenceIndex();
	}

private void throwError(final T loc) {
	String msg = "Illegal Order got " +loc.getContig()+":"+loc.getStart()+" after "+ this.prev_contig+":"+this.prev_start;
	if(this.strictDictOrder) {
		msg += ". Order must conform the sequence dictionary of the reference sequence.";
		}
	throw new IllegalArgumentException(msg);
	}

@Override
public T apply(final T t) {
	if(t==null) return null;
	final String contig = t.getContig();
	if (this.prev_contig==null) {
		if(this.dict!=null) {
			getTid(contig);
			}
		this.prev_contig = t.getContig();
		this.prev_start = t.getStart();
		if(!this.strictDictOrder) this.contig_seen.add(contig);
		}
	else if(!this.prev_contig.equals(contig)) {
		if (this.strictDictOrder ) {
			final int prev_tid = getTid(this.prev_contig);
			final int tid = getTid(contig);
			if (prev_tid > tid) {
				throwError(t);
				}
			}
		else {
			if(this.dict!=null) getTid(contig);
			if(this.contig_seen.contains(contig)) {
				throw new IllegalStateException("Bad order. Saw contig "+contig+" twice. Input must be sorted.");
				}
			this.contig_seen.add(contig);
			}
		this.prev_contig = t.getContig();
		this.prev_start = t.getStart();		
		}
	else if (this.prev_start > t.getStart()) {
		throwError(t);
		}
	else {
		this.prev_start = t.getStart();
		}
	return t;
	}

@Override
public OrderChecker<T> clone() {
	return new OrderChecker<>(this.dict, this.strictDictOrder);
	}

@Override
public String toString() {
	return "OrderChecker:"+ (this.dict==null?"without":"with")+" dict. strict order:"+this.strictDictOrder;
	}
}
