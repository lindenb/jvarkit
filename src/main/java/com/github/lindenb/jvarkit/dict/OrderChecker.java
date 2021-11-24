package com.github.lindenb.jvarkit.dict;

import java.util.HashSet;
import java.util.Set;
import java.util.function.UnaryOperator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;

public class OrderChecker<T extends Locatable> implements UnaryOperator<T>{
private final Set<String> contig_seen = new HashSet<>();
private final SAMSequenceDictionary dict;
private final Delegate delegate;
private abstract class Delegate implements UnaryOperator<T> {
	final Set<String> contig_seen = new HashSet<>();
	String prev_contig = null;
	int prev_start = -1;
	@Override
	public T apply(final T t) {
		if(t==null) return null;
		final String contig = t.getContig();
		if(prev_contig==null) {
			
			}
		return null;
	}
}

public OrderChecker(final SAMSequenceDictionary dict,boolean strictDictOrder) {
	this.dict = dict;
	}
	
@Override
public T apply(T t) {
	return this.delegate.apply(t);
	}
}
