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

/** utility to find a variant in a sam recor
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
		final StringBuilder sbr = new StringBuilder(ctx.getReference().getDisplayString());
		final StringBuilder sba = new StringBuilder(alt.getDisplayString());
		// REF ATTT
		// ALT AT
		while(	sba.length()>0 && 
				sbr.length()> sba.length() && 
				sbr.charAt(0)==sba.charAt(0)) {
				this.refpos++;
				sbr.deleteCharAt(0);
				sba.deleteCharAt(0);
			}
		this.del_size = sbr.length() - sba.length();
	}
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
	
	final List<Allele> ins = alts.stream().
			filter(A->A.length()> ref_allele_len).
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
					if(clen==norm.del_size && norm.refpos==ref1+1 && norm.del_size<0 /* TODO */) {
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
				for(final Allele a: ins) {
					final int ins_length = a.length() - ref_allele_len ;
					if(clen==ins_length && ctx.getStart()+ ref_allele_len==ref1) {
						final byte[] allele_bases = a.getBases();
						int j=0;
						for(j=0;j< ins_length;++j) {
							int b1x= read0+ref_allele_len+j;
							int b2x= ref_allele_len+j;
							if(b1x>=bases.length) break;
							if(b2x>=allele_bases.length) break;
							// read base
							final byte b1 = bases[b1x];
							// allele base
							final byte b2 = allele_bases[b2x];
							if (b1!=b2) {
								break;
								}
							}
						if(j==ins_length) {
							match.allele = Optional.of(a);
							match.read_pos = read0 ;
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
