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
package com.github.lindenb.jvarkit.samtools;

import java.util.ArrayList;
import java.util.List;
import java.util.function.UnaryOperator;

import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class SAMRecordLeftAligner implements UnaryOperator<SAMRecord> {
	public static final String OPT_DESC="Left-aligns any indels in the read data contained in a BAM or CRAM file. The same indel can often be placed at multiple positions and still represent the same haplotype. While it is a commonly used convention to place an indel at the left-most position";
	private GenomicSequence genomicSeq=null;
	private final ReferenceSequenceFile referenceSequenceFile;
	private boolean debug;
	private boolean previous_was_realigned=false;
	public SAMRecordLeftAligner(final ReferenceSequenceFile referenceSequenceFile) {
		this.referenceSequenceFile=referenceSequenceFile;
		}
	
	/** return true if previous record was realigned */
	public boolean previousRecordWasRealigned() {
		return this.previous_was_realigned;
		}
	
	private boolean isMatch(final CigarOperator op) {
		return op.equals(CigarOperator.M) || op.equals(CigarOperator.EQ);
		}
	
	private boolean isDel(final CigarOperator op) {
		return op.equals(CigarOperator.D) || op.equals(CigarOperator.N);
		}
	
	public void setDebug(boolean debug) {
		this.debug = debug;
		}
	
	public boolean isDebugging() {
		return debug;
		}
	
	@Override
	public SAMRecord apply(final SAMRecord rec) {
		previous_was_realigned=false;
		//if(!rec.getReadName().equals("RF02_1926_2453_0:0:2_1:0:0_e") || rec.getStart()!=1926)  return rec;
		if(rec.getReadUnmappedFlag() || this.referenceSequenceFile==null) return rec;
		final Cigar cigar = rec.getCigar();
		if(cigar==null || cigar.numCigarElements()<2 || cigar.getCigarElements().stream().noneMatch(OP->isDel(OP.getOperator()))) {
			return rec;
			}
		final byte[] bases = rec.getReadBases();
		if(bases==SAMRecord.NULL_SEQUENCE) return rec;
		boolean dirty=false;
		final List<CigarElement> L=new ArrayList<>(cigar.getCigarElements());
		int refpos1 = rec.getAlignmentStart();
		int readpos0  = 0;
		for(int i=0;i< L.size();++i) {
			final CigarElement ce = L.get(i);
			final  CigarOperator op = ce.getOperator();
			final int cigar_length = ce.getLength();
			
			switch(op) {
				case P: break;
				case H: break;
				case I: readpos0 += cigar_length; break;
				case S: readpos0 += cigar_length; break;
				case M: case X: case EQ:
					refpos1 += cigar_length;
					readpos0 += cigar_length;
					break;
				case D: case N:
					{
					int next_refpos1 = refpos1 + cigar_length;
					if(i>0 && i+1 < L.size()) {
						final CigarElement prev = L.get(i-1);
						final CigarElement next = L.get(i+1);
						final int next_length = next.getLength();
						final int next_readpos0 = readpos0;//no change
						if(isMatch(prev.getOperator()) && isMatch(next.getOperator()))
							{
							if(genomicSeq==null || !genomicSeq.hasName(rec.getReferenceName())) {
								genomicSeq = new GenomicSequence(this.referenceSequenceFile, rec.getReferenceName());
								}
							/* find bases in deletions where bases are the same than in ref */
							int extend=0;
							for(int x=0;
								x < cigar_length &&
								x< next_length -1 && /* -1 to avoid to have and empty CigarElement */
								next_readpos0 + x + 1 < bases.length ; /* we want one more base to test for mismatch after the correction see below*/
								x++) {
								final char base_in_del = Character.toUpperCase(genomicSeq.charAt(refpos1 + x- 1));
								final char base_in_next_read = (char)Character.toUpperCase(bases[next_readpos0 + x]);
								if(base_in_del!=base_in_next_read) break;
								if(isDebugging() ) {
									System.err.println(
											"extend\t"+
											rec.getPairedReadName()+
												"\tCONTIG:"+rec.getContig()+
												"\tREF["+(refpos1+x)+"]="+base_in_del+
												"\tBASE["+(next_readpos0 + x)+"]="+base_in_next_read+
												"\textend:"+extend
												);
									}
								extend++;
								}
							
							/* ok, to the other side, because we do not want the next Cigar element starting with a mismatch */
							while(extend>0) {
								final char base_in_del = Character.toUpperCase(genomicSeq.charAt((next_refpos1 + extend- 1) + 1));
								final char base_in_next_read = (char)Character.toUpperCase(bases[next_readpos0 + extend + 1]);
								/* it's ok, it's a match we can break here */
								if(base_in_del==base_in_next_read) break;
								if(isDebugging() ) {
									System.err.println(
											"shrink\t"+
											rec.getPairedReadName()+
												"\tCONTIG:"+rec.getContig()+
												"\tREF["+(refpos1+extend)+"]="+base_in_del+
												"\tBASE["+(next_readpos0 + extend)+"]="+base_in_next_read+
												"\textend:"+extend
												);
									}
								extend--;
								}
							
							if(extend>0) {
								dirty=true;
								L.set(i-1, new CigarElement(prev.getLength()+extend, prev.getOperator()));
								L.set(i+1, new CigarElement(next.getLength()-extend, next.getOperator()));
								next_refpos1+=extend;
								}
							
							}
						}
					refpos1 = next_refpos1;
					break;
					}
				}
			}
		
		if(!dirty) {
			return rec;
			}
		
		final Cigar newCigar = new Cigar(L);
		//System.err.println("Changed "+cigar.toString()+" to "+newCigar);
		previous_was_realigned=true;
		rec.setCigar(newCigar);
		return rec;
		}
	
	/** dispose GenomicSequence */
	public void dispose() {
		previous_was_realigned=false;
		this.genomicSeq=null;
		}
	}
