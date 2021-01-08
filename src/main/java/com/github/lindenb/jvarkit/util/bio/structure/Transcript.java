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
package com.github.lindenb.jvarkit.util.bio.structure;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.OptionalInt;

import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;

import htsjdk.samtools.util.Locatable;


public interface Transcript extends StrandedLocatable {
	/** get property map */
	public Map<String, String> getProperties();
	/** get associated gene */
	public Gene getGene();
	/** get transcript id */
	public String getId();
	/** return i-th exon scanning from 5' to 3', whatever strand */
	public Exon getExon(int index0);
	/** return i-th exon start scanning from 5' to 3', whatever strand */
	public int getExonStart(int index0);
	/** return i-th exon end scanning from 5' to 3', whatever strand */
	public int getExonEnd(int index0);
	/** return exons scanning from 5' to 3', whatever strand */
	public List<Exon> getExons();
    public char getStrand();
    public int getExonCount();
    /** weird transcripts when no exons was defined */
    public default boolean hasExon() {
    	return getExonCount()>0;
    }
    /** get all CDS in this exon */
	public List<Cds> getAllCds();

	/** returns true hasCodonStrop && hasCodonStart */
    public default boolean hasCDS() {
    	return hasCodonStartDefined() && hasCodonStopDefined();
    }
	
	/** returns true if getIntronCount()>0 */
    public default boolean hasIntron() {
    	return getIntronCount()>0;
    }
    /** return the number of introns */
    public int getIntronCount();
	public List<Intron> getIntrons();
	/** return i-th intron scanning from 5' to 3', whatever strand */
	public Intron getIntron(int index0);

	/** get genomic transcription start position in the genome */
    public int getTxStart();
	/** get genomic transcription end position in the genome */
    public int getTxEnd();

    public boolean isNonCoding();
	public boolean isCoding();
	
	
	/** get transcript length (cumulative exons sizes )*/
	public default int getTranscriptLength() {
		return getExons().stream().
				mapToInt(T->T.getLengthOnReference()).
				sum();
		}
	
	/** get transcript length (cumulative CDS sizes )*/
	public default int getCodingDNALength() {
		return getAllCds().stream().
				mapToInt(T->T.getLengthOnReference()).
				sum();
		}
	
	/** return UTR on 5' side of genomic reference */
	public Optional<UTR> getUTR5();
	/** return UTR on 3' side of genomic reference */
	public Optional<UTR> getUTR3();

	/** return UTR on 5' side of transcript. strand matters */
	public default Optional<UTR> getTranscriptUTR5() {
		if(!hasCodonStartDefined()) return Optional.empty();
		if(isPositiveStrand()) return getUTR5();
		if(isNegativeStrand()) return getUTR3();
		return Optional.empty();
	}
	
	/** return UTR on 5' side of transcript. strand matters */
	public default Optional<UTR> getTranscriptUTR3() {
		if(!hasCodonStopDefined()) return Optional.empty();
		if(isPositiveStrand()) return getUTR3();
		if(isNegativeStrand()) return getUTR5();
		return Optional.empty();
	}

	/** return a list of UTRs, may be empty */
	public default List<UTR> getUTRs() {
		final Optional<UTR> utr5 =  getUTR5();
		final Optional<UTR> utr3 =  getUTR3();
		if(!utr5.isPresent() && !utr3.isPresent()) return Collections.emptyList();
		if(utr5.isPresent() && !utr3.isPresent()) return Collections.singletonList(utr5.get());
		if(!utr5.isPresent() && utr3.isPresent()) return Collections.singletonList(utr3.get());
		return Arrays.asList(utr5.get(),utr3.get());
	}
	
	
	
	public default boolean hasUTR() {
		return getUTR5().isPresent()|| getUTR3().isPresent();
	}
	
	/** return start-codon if (+) strand . return stop-codon if (-) strand */
	public default Optional<Codon> getLeftmostCodon() {
		if(isPositiveStrand()) return getCodonStart();
		if(isNegativeStrand()) return getCodonStop();
		return Optional.empty();
	}
	
	/** return stop-codon if (+) strand . return start-codon if (-) strand */
	public default Optional<Codon> getRightmostCodon() {
		if(isPositiveStrand()) return getCodonStop();
		if(isNegativeStrand()) return getCodonStart();
		return Optional.empty();
	}
	
	/** get Codon start */
	public Optional<Codon> getCodonStart();
	/** get Codon end */
	public Optional<Codon> getCodonStop();
	
	/** test if start_codon defined. eg:  wget -q -O - "http://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz" | gunzip -c | grep ENSG00000060138 */
	public boolean hasCodonStartDefined();	
	/** test if stop codon defined. eg:  wget -q -O - "http://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz" | gunzip -c | grep ENST00000327956 | grep stop */
	public boolean hasCodonStopDefined();

	/**return genomic interval spanning between codon start and codon stop if they both exists*/
	public default Optional<Locatable> getCdsInterval() {
		if(!hasCodonStartDefined() || !hasCodonStopDefined()) return Optional.empty();
		int x1=this.getEnd();
		int x2=this.getStart();

		Locatable codon = getCodonStart().get();
		x1=Math.min(codon.getStart(), x1);
		x2=Math.max(codon.getEnd(), x2);
		
		codon = getCodonStop().get();
		x1=Math.min(codon.getStart(), x1);
		x2=Math.max(codon.getEnd(), x2);
		return Optional.of(new SimpleInterval(getContig(), x1, x2));
		}
	
	/** 1-based leftmost position of the cds in genomic order (equivalent of ucsc knownGene.getCdsStart) */
	public default OptionalInt getCdsStart() {
		final Integer p= getCdsInterval().map(A->A.getStart()).orElse(null);
		return p==null?OptionalInt.empty() : OptionalInt.of(p.intValue());
		}
	
	/** 1-based rightmost position of the cds in genomic order (equivalent of ucsc knownGene.getCdsStart) */
	public default OptionalInt getCdsEnd() {
		final Integer p= getCdsInterval().map(A->A.getEnd()).orElse(null);
		return p==null?OptionalInt.empty() : OptionalInt.of(p.intValue());
		}

}
