/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.function.BiFunction;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

/** utility to find a variant in a sam record
 * works best for indel if VCF was normalized
 * 
 * */
public class FindVariantInSamRecord implements BiFunction<SAMRecord, VariantContext, FindVariantInSamRecord.Match> {
	private boolean use_clip =false;

public static interface Match {
	public SAMRecord getRecord();
	public VariantContext getVariant();
	public Optional<Allele> getAllele();
	}
private static class MatchImpl implements Match {
	private final SAMRecord rec;
	private final VariantContext ctx;
	private Optional<Allele> allele = Optional.empty();
	private int read_pos = -1;
	MatchImpl(final SAMRecord rec,final VariantContext ctx) {
		this.rec = rec;
		this.ctx = ctx;
		}
	@Override
	public SAMRecord getRecord() {
		return rec;
		}
	@Override
	public VariantContext getVariant() {
		return ctx;
		}
	public Optional<Allele> getAllele() {
		return allele;
		}
	}

public FindVariantInSamRecord setUseClip(boolean use_clip) {
	this.use_clip = use_clip;
	return this;
	}

public boolean isUseClip() {
	return use_clip;
	}

@Override
public final Match apply(final SAMRecord t,final VariantContext u) {
	return find(t,u);
	}

private boolean isValidAllele(final Allele a) {
	return AcidNucleics.isATGC(a);
}

private static class NormalizedDelVariant {
	final Allele alt;
	int refpos;
	final int del_size;
	NormalizedDelVariant(final VariantContext ctx,final Allele alt) {
		this.alt = alt;
		this.refpos = ctx.getStart();
		final byte[] sr = ctx.getReference().getBases();
		final byte[] sa = alt.getBases();
		int x = 0;
		// REF ATTT
		// ALT A
		while(sr.length > sa.length &&
			x < sa.length &&
			sr[x]==sa[x]) {
			this.refpos++;
			x++;
			}
		this.del_size = sr.length - x;
	}
}


private static class NormalizedInsVariant {
	final Allele alt;
	int refpos;
	final int dx;
	final byte[] sa;
	NormalizedInsVariant(final VariantContext ctx,final Allele alt) {
		this.alt = alt;
		this.refpos = ctx.getStart();
		final byte[] sr = ctx.getReference().getBases();
		this.sa = alt.getBases();
		if(this.sa.length <= sr.length) throw new IllegalArgumentException();
		int x = 0;
		// REF A
		// ALT ATTT
		while(sr.length < this.sa.length &&
			x < sr.length &&
			sr[x]==this.sa[x]) {
			this.refpos++;
			x++;
			}
		this.dx=x;
	}
	int size() { return this.sa.length - this.dx;}
	byte at(int idx) { return this.sa[idx+this.dx];}
}




public Match find(final SAMRecord record,final VariantContext ctx) {
	final MatchImpl match = new MatchImpl(record,ctx);

	final Locatable recloc = isUseClip()?
			new SimpleInterval(record.getContig(),record.getUnclippedStart(),record.getUnclippedEnd()):
			record
			;
	if(!ctx.overlaps(recloc)) {
		return match;
	}
	final byte[] bases = record.getReadBases();
	if(bases==null || bases==SAMRecord.NULL_SEQUENCE) {
		return match;
		}
	final Cigar cigar = record.getCigar();
	if(cigar==null || cigar.isEmpty()) {
		return match;
		}
	
	if(!isValidAllele(ctx.getReference())) {
		return match;
	}
	
	final List<Allele> alts = ctx.getAlternateAlleles().
			stream().
			filter(A->isValidAllele(A)).
			collect(Collectors.toList());
	
	if(alts.isEmpty()) {
		return match;
	}
	
	final int ref_allele_len = ctx.getReference().length();
	
	final List<NormalizedDelVariant> dels = alts.stream().
			filter(A->A.length()< ref_allele_len).
			map(A->new NormalizedDelVariant(ctx,A)).
			collect(Collectors.toList());
	
	final List<NormalizedInsVariant> ins = alts.stream().
			filter(A->A.length()> ref_allele_len).
			map(A->new NormalizedInsVariant(ctx,A)).
			collect(Collectors.toList());

	
	
	final List<Allele> subst = new ArrayList<>(ctx.getNAlleles());
	//add ref
	subst.add(ctx.getReference());
	// add alleles having same length as ref
	subst.addAll(alts.stream().
			filter(A->A.length()==ref_allele_len).
			collect(Collectors.toList())
			);

		
	int read0= 0 ;
	int ref1 = record.getUnclippedStart();
	for(final CigarElement ce:cigar) {
		if(ref1 > ctx.getEnd()) break;
		final CigarOperator op = ce.getOperator();
		final int clen = ce.getLength();
		
		switch(op) {
			case P:break;
			case H: ref1 += clen;break;
			case D: case N: {
				// REF AAAAAAAAAAAAAAA
				// ALT AAAA<-- len -->
				for(final NormalizedDelVariant norm: dels) {
					//System.err.println("ref1:"+ref1+"=="+norm.refpos+" ctx.start="+ctx.getContig()+":"+ctx.getStart()+":"+ctx.getReference()+" norm.pos="+norm.refpos+" norm.size:"+norm.del_size+"=="+clen+" alt:"+norm.alt+" "+record.getStart()+":"+record.getCigarString());
					if(clen==norm.del_size && norm.refpos==ref1) {
						//System.err.println("OK");
						match.allele = Optional.of(norm.alt);
						match.read_pos = read0 ;
						return match;
						}
					}
				ref1 += clen;
				break;
				}
			case I: {
				// REF AAAA<-- len -->
				// ALT AAAAAAAAAAAAAAA
				for(final NormalizedInsVariant norm: ins) {
					//System.err.println("INS ref1:"+ref1+"=="+norm.refpos+" ctx.start="+ctx.getContig()+":"+ctx.getStart()+":"+ctx.getReference()+" norm.pos="+norm.refpos+" norm.size:"+norm.size()+"=="+clen+" alt:"+norm.alt+" "+record.getStart()+":"+record.getCigarString());
					if(clen==norm.size() && norm.refpos==ref1) {
						int j=0;
						for(j=0;j< clen;++j) {
							int b1x= read0 +j;
							if(b1x>= bases.length) break;
							// read base
							final byte b1 = bases[b1x];
							// allele base
							final byte b2 = norm.at(j);
							if (b1!=b2) {
								break;
								}
							}
						if(j==clen) {
							match.allele = Optional.of(norm.alt);
							return match;
							}
						}
					}
				read0 += clen;
				break;
				}
			case M: case X: case EQ: case S: {
				if(op.equals(CigarOperator.S) && !isUseClip()) {
					read0 += clen;
					ref1 += clen;
					break;
					}
				for(int x=0;x + ref_allele_len <= clen;++x) {
					if(ref1+x!=ctx.getStart()) continue;
					for(final Allele a: subst) {

						final byte[] allele_bases = a.getBases();
						if(allele_bases.length!=ref_allele_len) throw new IllegalStateException();
						int j=0;
						for(j=0;j< ref_allele_len ;++j) {
							// read base
							final byte b1 = bases[read0+x+j];
							// allele base
							final byte b2 = allele_bases[j];
							if (b1!=b2) {
								break;
								}
							}
						if(j== ref_allele_len) {
							match.allele = Optional.of(a);
							match.read_pos = read0+x;
							return match;
							}
						}
					}
				read0 += clen;
				ref1 += clen;
				break;
				}
			default: throw new IllegalArgumentException(op.toString());
			}
		}
	
	return match;
	}
}
