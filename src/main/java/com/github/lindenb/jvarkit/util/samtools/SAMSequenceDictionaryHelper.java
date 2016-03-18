/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
* 2015 creation

*/
package com.github.lindenb.jvarkit.util.samtools;

import java.util.Optional;
import java.util.stream.Collectors;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.vcf.VCFHeader;

public class SAMSequenceDictionaryHelper {
private final SAMSequenceDictionary dict;
private boolean throwingOnError = true;
private boolean trimName = false;
private boolean searchAltName = false;

public SAMSequenceDictionaryHelper(final SAMSequenceDictionary dict) {
	this.dict=dict;
}



public SAMSequenceDictionaryHelper(final SAMSequenceDictionaryHelper cp) {
	this(cp.dict);
	this.throwingOnError = cp.throwingOnError;
	this.trimName = cp.trimName;
	this.searchAltName = cp.searchAltName;
}

public SAMSequenceDictionaryHelper(final VCFHeader header) {
	if(header==null) throw new IllegalArgumentException("VCF header is null");
	if(header.getSequenceDictionary() == null ) throw new IllegalArgumentException("VCF header is missing a VCF dictionary ( in the header , lines starting with ##contig)");
	this.dict = header.getSequenceDictionary();
}

public SAMSequenceDictionaryHelper(final SAMFileHeader header) {
	if(header==null) throw new IllegalArgumentException("SAM File header is null");
	if(header.getSequenceDictionary() == null ) throw new IllegalArgumentException("VCF header is missing a VCF dictionary ( in the header , lines starting with @SEQ)");
	this.dict = header.getSequenceDictionary();
}

public SAMSequenceDictionaryHelper() {
	this(new SAMSequenceDictionary());
}

public SAMSequenceDictionary getDictionary() {
	return dict;
}

public double fraction(final String contig,int pos) {
	return scannedSoFar(contig,pos)/(double)getDictionary().getReferenceLength();
}


public long scannedSoFar(final String contig,int pos) {
	long sum=0L;
	for(final SAMSequenceRecord sr:getDictionary().getSequences()) {
		if(sr.getSequenceName().equals(contig)) return sum+pos;
		sum += sr.getSequenceLength();
	}
	return 0L;
}

private String trim(final String s) {
	return s==null || !this.trimName?s:s.trim();
}

public boolean containsContig(final String userContig) {
	if(userContig==null || userContig.isEmpty()) return false;
	return findSAMSequenceRecord(false,userContig).isPresent();
	
}

public void setSearchingAltName(boolean searchAltName) {
	this.searchAltName = searchAltName;
}

public boolean isSearchingAltName() {
	return searchAltName;
}

public Optional<SAMSequenceRecord> findSAMSequenceRecord(final String userContig) { 
	return findSAMSequenceRecord(isThrowingOnError(),userContig);
}

private Optional<SAMSequenceRecord> findSAMSequenceRecord(boolean exception,final String userContig) {
	String contig = trim(userContig);
	if(contig==null) return empty(exception,"Input string is null.");
	if(contig.isEmpty()) return empty(exception,"Input string is empty.");

	SAMSequenceRecord rec = getDictionary().getSequence(contig);
	if(rec!=null) return Optional.of(rec);

	if(isSearchingAltName()) {
		if(contig.startsWith("chr")) {
			rec = getDictionary().getSequence(contig.substring(3));
			if(rec!=null) return Optional.of(rec);
			if(contig.equals("chrM")) {
				rec = getDictionary().getSequence("MT");
				if(rec!=null) return Optional.of(rec);
			}
		} else {
			rec = getDictionary().getSequence("chr"+contig);
			if(rec!=null) return Optional.of(rec);
			if(contig.equals("MT")) {
				rec = getDictionary().getSequence("chrM");
				if(rec!=null) return Optional.of(rec);
			}
		}
	}
	if(!exception) return Optional.empty();//prevent building large string below
	return empty(exception,"Cannot find \""+contig+" \" in sequence dictionary. Available contigs are: "+
			getDictionary().getSequences().stream().map(C->C.getSequenceName()).collect(Collectors.joining(" ; "))
			);
	}

public Interval convert(final Interval interval) {
	if(interval==null)throw new IllegalArgumentException("Input interval is null.");
	Optional<SAMSequenceRecord> ssr = this.findSAMSequenceRecord(interval.getContig());
	if(!ssr.isPresent()) return null;
	return new Interval(
			ssr.get().getSequenceName(),
			interval.getStart(),
			interval.getEnd(),
			interval.isNegativeStrand(),
			interval.getName()
			);
	}

public Interval convert(final String contig,int start,int end) {
	return convert(new Interval(contig, start, end));
	}


public void setThrowingOnError(boolean throwingOnError) {
	this.throwingOnError = throwingOnError;
}

public boolean isThrowingOnError() {
	return throwingOnError;
}

private <T> Optional<T> empty(boolean exception,String message) {
	return empty(exception,message,null);
}

private <T> Optional<T> empty(boolean exception,String message,Throwable cause) {
	if(exception) {
		throw new IllegalArgumentException(message, cause);
	}
	return Optional.empty();
}



}
